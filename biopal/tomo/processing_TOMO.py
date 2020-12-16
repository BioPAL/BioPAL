import os
import logging
import numpy as np
from scipy.signal import medfilt2d
from gdalconst import GA_ReadOnly
from osgeo import gdal
from matplotlib import pylab as plt
from biopal.statistics.utility_statistics import (
    build_filtering_matrix,
    Covariance2D2Correlation2D,
    SKPD_processing,
    main_correlation_estimation_SR,
    covariance_matrix_vec2mat,
)


def UpperThresholdForestHeight(power_cube, opt_str):
    """It finds the last drop of the power spectrum from opt_str.thr times the
    % peak of the spectrum. The corresponding elevation is returned.
    %
    % INPUT
    %      power_cube: [Nr x Nc x Nz] positive real valued cube of
    %                  backscattered power
    %      opt_str.
    %              thr: threshold for the spectrum drop (between 0 and 1)
    %              z: [Nz x 1] elevation axis
    % OUTPUT
    %       out_str.
    %               z: [Nr x Nc] canopy elevation map
    %
    % DEPENDENCIES
    %             FindStep3"""
    Nr, Na, Nz = power_cube.shape
    power_peak = np.max(power_cube, axis=2)
    greater_power_mask = power_cube > power_peak.reshape((Nr, Na, 1)) * opt_str.thr

    last_drop_linear_indices = FindStep3(greater_power_mask, 'last')

    class out_str:
        pass

    out_str.z = opt_str.z[last_drop_linear_indices]
    out_str.peak = power_peak
    return out_str


def FindStep3(I, varargin):
    """This routin finds the first or the last step in a logical 3D along the
    % third direction. The first step is a rising one whereas the last is a
    % drop.
    %
    % INPUT
    %      I: [Nr x Nc x N] logical array to be explored
    %      varargin{1}: (optional. Default: 'first') either 'first' or 'last'
    %
    % OUTPUT
    %       linear_indices: [Nr x Nc] matrix of linear indices where the step
    %                       has been detected"""
    Nr, Na, Nz = I.shape
    I = I > 0

    if varargin == 'first':

        rise_mask = 1 - (np.roll(I, (0, 0, 1), (0, 1, 2))) & I
        rise_mask[:, :, 0] = False
        rise_mask[:, :, -1] = False

        rise_mask_cumsum = np.cumsum(rise_mask, axis=2)

        linear_indices = np.sum(
            ((rise_mask_cumsum == 1) & rise_mask) * np.arange(I.shape[2]).reshape(1, 1, I.shape[2]),
            axis=2,
        )

    elif varargin == 'last':

        drop_mask = 1 - (np.roll(I, (0, 0, -1), (0, 1, 2))) & I
        drop_mask[:, :, 0] = False
        drop_mask[:, :, -1] = False

        drop_mask_cumsum = np.cumsum(drop_mask, axis=2)

        linear_indices = np.sum(
            ((drop_mask_cumsum == drop_mask_cumsum[:, :, -1].reshape((Nr, Na, 1))) & drop_mask)
            * np.arange(I.shape[2]).reshape(1, 1, I.shape[2]),
            axis=2,
        )

    else:

        print('Unrecognized input, aborting.\n')
        linear_indices = []

    return linear_indices


def BiomassForestHeightSKPD(
    data_stack,
    cov_est_window_size,
    pixel_spacing_slant_rg,
    pixel_spacing_az,
    incidence_angle_rad,
    kz_stack,
    vertical_vector,
    proc_conf,
):

    power_threshold = proc_conf.power_threshold

    # data_stack is a dictionary of two nested dictionaries composed as:
    # data_stack[ acquisition_name ][ polarization ]

    num_acq = len(data_stack)
    acq_names = list(data_stack.keys())
    first_acq_dict = data_stack[acq_names[0]]
    pol_names = list(first_acq_dict.keys())
    num_pols = len(pol_names)
    Nrg, Naz = first_acq_dict[pol_names[0]].shape
    Nz = np.size(vertical_vector)

    # Covariance estimation
    MPMB_correlation, rg_vec_subs, az_vec_subs, subs_F_r, subs_F_a = main_correlation_estimation_SR(
        data_stack,
        cov_est_window_size,
        pixel_spacing_slant_rg,
        pixel_spacing_az,
        incidence_angle_rad,
    )

    Nrg_subs = rg_vec_subs.size
    Naz_subs = az_vec_subs.size

    # Initialization of the SKPD routine
    class SKPD_kernel_opt_str:
        pass

    SKPD_kernel_opt_str.num_acq = num_acq
    SKPD_kernel_opt_str.num_pols = num_pols
    SKPD_kernel_opt_str.Nsubspaces = 2
    SKPD_kernel_opt_str.Nparam = 2
    SKPD_kernel_opt_str.error = np.zeros((Nrg_subs, Naz_subs, 4, 2))

    # Single polarimetric channel selector
    wpol = (np.kron(np.eye(num_pols), np.ones((1, num_acq)))) > 0

    tomo_cube = np.zeros((Nrg_subs, Naz_subs, Nz))
    Nrg_subs_string = str(Nrg_subs)
    for rg_sub_idx in np.arange(Nrg_subs):

        logging.info('   Heigth step ' + str(rg_sub_idx + 1) + ' of ' + Nrg_subs_string)
        for az_sub_idx in np.arange(Naz_subs):

            # Spectra estimation initialization
            class spectra:
                pass

            spectra.temp = np.zeros((Nz, 4))
            current_MPMB_correlation = MPMB_correlation[:, :, rg_sub_idx, az_sub_idx]

            # SKPD processing
            SKPD_kernel_out_str = SKPD_processing(current_MPMB_correlation, SKPD_kernel_opt_str)

            if np.any(SKPD_kernel_out_str.error):
                # Calibrating with respect to the linked phases of the
                # best polarimetric channel
                Rcoh_thin = Covariance2D2Correlation2D(
                    np.reshape(
                        current_MPMB_correlation[
                            wpol[0, :][:, np.newaxis] * wpol[0, :][np.newaxis, :]
                        ],
                        (num_acq, num_acq),
                    )
                )
                Rcoh_fat = Covariance2D2Correlation2D(
                    np.reshape(
                        current_MPMB_correlation[
                            wpol[1, :][:, np.newaxis] * wpol[1, :][np.newaxis, :]
                        ],
                        (num_acq, num_acq),
                    )
                )
                Rincoh_thin = Covariance2D2Correlation2D(
                    np.reshape(
                        current_MPMB_correlation[
                            wpol[2, :][:, np.newaxis] * wpol[2, :][np.newaxis, :]
                        ],
                        (num_acq, num_acq),
                    )
                )
                Rincoh_fat = np.random.randn(num_acq, num_acq) + 1j * np.random.randn(
                    num_acq, num_acq
                )

            else:
                # Scattering mechanisms
                Rcoh_thin = Covariance2D2Correlation2D(SKPD_kernel_out_str.Rcoh_thin)
                Rcoh_fat = Covariance2D2Correlation2D(SKPD_kernel_out_str.Rcoh_fat)
                Rincoh_thin = Covariance2D2Correlation2D(SKPD_kernel_out_str.Rincoh_thin)
                Rincoh_fat = Covariance2D2Correlation2D(SKPD_kernel_out_str.Rincoh_fat)

            # Steering matrix
            current_kz = np.zeros((num_acq, 1))
            for b_idx, stack_curr in enumerate(kz_stack.values()):
                current_kz[b_idx] = stack_curr[rg_vec_subs[rg_sub_idx], az_vec_subs[az_sub_idx]]

            A = np.exp(1j * current_kz * vertical_vector) / num_acq

            # Spectra estimation
            for m in np.arange(4):
                currR = (
                    (m == 0) * Rcoh_thin
                    + (m == 1) * Rcoh_fat
                    + (m == 2) * Rincoh_thin
                    + (m == 3) * Rincoh_fat
                )
                if proc_conf.enable_super_resolution:
                    # Capon
                    currR = currR + proc_conf.regularization_noise_factor * np.eye(currR.shape[0])
                    spectra.temp[:, m] = 1 / np.diag(
                        np.abs(A.conj().transpose() @ np.linalg.inv(currR) @ A)
                    )
                else:
                    spectra.temp[:, m] = np.diag(np.abs(A.conj().transpose() @ currR @ A))

            # Volume mechanism recognized thanks to its higher elevation
            max_index = np.argmax(spectra.temp, axis=0)
            max_m = np.argmax(max_index)
            tomo_cube[rg_sub_idx, az_sub_idx, :] = spectra.temp[:, max_m]

    # Estimating canopy elevation by looking at the decay
    class opt_str:
        pass

    opt_str.z = vertical_vector
    opt_str.thr = power_threshold  # Power decay threshold With respect to the peak value
    out_str = UpperThresholdForestHeight(tomo_cube, opt_str)
    canopy_height = out_str.z
    power_peak = out_str.peak
    canopy_height = medfilt2d(canopy_height.astype('float64'), kernel_size=5)

    return canopy_height, power_peak, rg_vec_subs, az_vec_subs, subs_F_r, subs_F_a


def tomo_cube_generator(z, I, model_sign, kz_stack):
    # Beamforming (Fourier spectrum estimation)
    Nz = np.size(z)
    Nr, Na, N = I.shape
    # Initialization
    tomo_cube = np.zeros((Nr, Na, Nz), dtype=np.complex128)

    print('        Beamforming...')
    for z_ind in np.arange(Nz):
        steering_array = np.exp(-1j * model_sign * kz_stack * z[z_ind]) / np.sqrt(N)
        tomo_cube[:, :, z_ind] = np.sum(I * steering_array, axis=2)
    print('        done.')
    return tomo_cube


def tomo_power_computation_OLD(savPath, win_x, win_y, pol_string, save_flag=False):
    loadPath = os.path.abspath(savPath)
    loadlist = os.listdir(loadPath)
    for ind in loadlist:
        if 'x_ax' in ind:
            x_ax = np.load(os.path.join(loadPath, ind))
            Na = np.size(x_ax)
            da = x_ax[1] - x_ax[0]
        if 'y_ax' in ind:
            y_ax = np.load(os.path.join(loadPath, ind))
            Nr = np.size(y_ax)
            dr = y_ax[1] - y_ax[0]
        if 'z_ax' in ind:
            z_ax = np.load(os.path.join(loadPath, ind))
            Nz = np.size(z_ax)

    # Range filter matrix
    class filtering_matrix_opt_str:
        pass

    NW = np.round(win_y / dr).astype('int64')
    filtering_matrix_opt_str.Nin = Nr
    filtering_matrix_opt_str.NWmono = np.floor((NW - 1) / 2).astype('int64')
    filtering_matrix_opt_str.subs_F = filtering_matrix_opt_str.NWmono
    Fr, rg_subs_ind, Rnorm = build_filtering_matrix(filtering_matrix_opt_str)
    Nrg_subs = rg_subs_ind.size
    Fr_normalized = Rnorm @ Fr

    # Azimuth filter matrix
    NW = np.round(win_x / da).astype('int64')
    filtering_matrix_opt_str.Nin = Na
    filtering_matrix_opt_str.NWmono = np.floor((NW - 1) / 2).astype('int64')
    filtering_matrix_opt_str.subs_F = filtering_matrix_opt_str.NWmono
    Fa, az_subs_ind, Anorm = build_filtering_matrix(filtering_matrix_opt_str)
    Naz_subs = az_subs_ind.size
    Fa_normalized = (Anorm @ Fa).T

    # Initialization
    power_cube = np.zeros((Nrg_subs, Naz_subs, Nz))
    n = 0
    ## Backscattered power estimation
    print('Power estimation...')
    counter = 0
    for fileName in loadlist:
        if '_tomo_' in fileName:
            fileName_split = fileName.split('_')
            current_pol_string = fileName_split[-3]
            slice_idx = int(fileName_split[-1].split('.')[0])

            if pol_string == current_pol_string:

                slice_loaded = np.load(os.path.join(loadPath, fileName))

                counter = counter + 1
                print('    adding slice ' + str(counter) + ' of ' + str(Nz))
                power_cube[:, :, slice_idx] = (
                    Fr_normalized @ np.abs(slice_loaded) ** 2 @ Fa_normalized
                )
                n += 1
                print('    {:d}'.format(int(n / Nz * 100)) + '%')

    x_ax = x_ax[az_subs_ind]
    y_ax = y_ax[rg_subs_ind]

    if save_flag:
        np.save(os.path.join(savPath, 'power_cube'), power_cube)
        np.save(os.path.join(savPath, 'x_ax_power_cube'), x_ax)
        np.save(os.path.join(savPath, 'x_ax_power_cube'), y_ax)
        np.save(os.path.join(savPath, 'x_ax_power_cube'), z_ax)

    return x_ax, y_ax, z_ax, power_cube


def tomo_power_computation(
    cube_sllices_folder, equi7_zone_name, equi7_tile_name, win_x, win_y, save_path=''
):

    start_a = 4460
    start_r = 3880
    Na = 140
    Nr = 160

    slices_folders = os.listdir(cube_sllices_folder)
    Nz = len(slices_folders)

    ## Backscattered power estimation
    print('Power estimation...')
    for slice_idx, slice_folder in enumerate(slices_folders):

        cube_slice_fname = os.path.join(
            cube_sllices_folder,
            slice_folder,
            equi7_zone_name,
            equi7_tile_name,
            'tomo_cube_slice_' + equi7_zone_name[6:] + '_' + equi7_tile_name + '.tif',
        )

        data_driver = gdal.Open(cube_slice_fname, GA_ReadOnly)
        for cov_idx in np.arange(6):
            if cov_idx == 0:
                cov_vec_curr = (
                    data_driver.GetRasterBand(int(cov_idx + 1))
                    .ReadAsArray(start_r, start_a, Nr, Na)
                    .T
                )
                plt.figure()
                plt.imshow(cov_vec_curr, cmap='jet')
                plt.pause(1)

                Nr = cov_vec_curr.shape[0]
                Na = cov_vec_curr.shape[1]
                cov_vec = np.zeros((6, Nr, Na))
                cov_vec[cov_idx, :, :] = cov_vec_curr
            else:
                cov_vec[cov_idx, :, :] = (
                    data_driver.GetRasterBand(int(cov_idx + 1))
                    .ReadAsArray(start_r, start_a, Nr, Na)
                    .T
                )

        print('reshape input cov...')
        cov_mat = covariance_matrix_vec2mat(cov_vec)
        print('done')

        if slice_idx == 0:
            # Initialization

            geotransform_equi7 = data_driver.GetGeoTransform()
            dr = np.abs(geotransform_equi7[1])
            da = np.abs(geotransform_equi7[5])

            # Range filter matrix
            class filtering_matrix_opt_str:
                pass

            NW = np.round(win_y / dr).astype('int64')
            filtering_matrix_opt_str.Nin = Nr
            filtering_matrix_opt_str.NWmono = np.floor((NW - 1) / 2).astype('int64')
            filtering_matrix_opt_str.subs_F = filtering_matrix_opt_str.NWmono
            Fr, rg_subs_ind, Rnorm = build_filtering_matrix(filtering_matrix_opt_str)
            Nrg_subs = rg_subs_ind.size
            Fr_normalized = Rnorm @ Fr

            # Azimuth filter matrix
            NW = np.round(win_x / da).astype('int64')
            filtering_matrix_opt_str.Nin = Na
            filtering_matrix_opt_str.NWmono = np.floor((NW - 1) / 2).astype('int64')
            filtering_matrix_opt_str.subs_F = filtering_matrix_opt_str.NWmono
            Fa, az_subs_ind, Anorm = build_filtering_matrix(filtering_matrix_opt_str)
            Naz_subs = az_subs_ind.size
            Fa_normalized = (Anorm @ Fa).T

            power_cube = np.zeros((Nrg_subs, Naz_subs, Nz))

        print('    adding ' + slice_folder)
        power_cube[:, :, slice_idx] = Fr_normalized @ np.abs(cov_mat) ** 2 @ Fa_normalized

    if save_path:
        np.save(os.path.join(save_path, 'power_cube'), power_cube)

    return power_cube


def tomo_plot(
    cube_slices_folder,
    pol_name,
    east_slice_index,
    north_slice_index,
    z_slice_index,
    equi7_zone_name,
    equi7_tile_name,
    plot_flag,
    pixel_extent='',
    save_power_cube_path='',
):

    if pol_name == 'hh':
        pol_idx = 0
    elif pol_name == 'hv':
        pol_idx = 3
    elif pol_name == 'vv':
        pol_idx = 5

    z_file = os.path.join(cube_slices_folder, 'vertical_vector.txt')
    if not os.path.exists(z_file):
        raise ValueError(
            'A file named "vertical_vector.txt" should be present in the slices folder '
            + cube_slices_folder
        )
    z_file_id = open(z_file, 'r')
    slices_folders = os.listdir(cube_slices_folder)
    num_z = len(slices_folders)
    z_ax = np.zeros(num_z)
    for slice_idx, cube_slice_curr in enumerate(slices_folders):
        if not 'tomo_cube_slice_' in cube_slice_curr:
            continue

        print('Reading slice {} of {}'.format(slice_idx + 1, num_z))
        line_curr = z_file_id.readline()
        idx_altitude = line_curr.find('altitude:')

        z_ax[slice_idx] = line_curr[idx_altitude + 10 : -5]

        cube_slice_fname = os.path.join(
            cube_slices_folder,
            cube_slice_curr,
            equi7_zone_name,
            equi7_tile_name,
            'tomo_cube_slice_' + equi7_zone_name[6:] + '_' + equi7_tile_name + '.tif',
        )
        print('cube_slice_fname: ', cube_slice_fname)
        data_driver = gdal.Open(cube_slice_fname, GA_ReadOnly)
        if slice_idx == 0:
            TEMP = data_driver.GetRasterBand(int(pol_idx + 1)).ReadAsArray()

            tomo_cube = np.zeros((num_z, TEMP.shape[0], TEMP.shape[1]))

            geotransform_equi7 = data_driver.GetGeoTransform()
            east_len = TEMP.shape[0]
            north_len = TEMP.shape[1]

            east_0 = geotransform_equi7[0]
            east_end = east_0 + geotransform_equi7[1] * east_len
            north_0 = geotransform_equi7[3]
            north_end = north_0 + geotransform_equi7[5] * north_len

            north_ax = np.linspace(0, north_len, num=north_len) * geotransform_equi7[5] + north_0
            east_ax = np.linspace(0, east_len, num=east_len) * geotransform_equi7[1] + east_0

        tomo_cube[slice_idx, :, :] = np.abs(
            data_driver.GetRasterBand(int(pol_idx + 1)).ReadAsArray()
        )

    if 1:
        if plot_flag == 'db':

            def fun_plot(temp_plot):
                return np.log10(temp_plot)

        else:

            def fun_plot(temp_plot):
                return temp_plot

        if pixel_extent:
            east_0 = east_ax[pixel_extent[0]] / 1000
            east_end = east_ax[pixel_extent[1]] / 1000
            north_0 = north_ax[pixel_extent[2]] / 1000
            north_end = north_ax[pixel_extent[3]] / 1000

        plt.figure()
        tomo_plt = np.squeeze(tomo_cube[z_slice_index, :, :])
        if pixel_extent:
            tomo_plt = tomo_plt[
                pixel_extent[0] : pixel_extent[1], pixel_extent[2] : pixel_extent[3]
            ]
        ax = plt.subplot(3, 1, 1)

        plt.imshow(
            fun_plot(tomo_plt.T),
            interpolation='none',
            origin='upper',
            extent=[east_0, east_end, north_end, north_0],
            cmap='jet',
        )
        plt.title('z = {:3.0f} [m]'.format(z_ax[z_slice_index]))
        plt.xlabel('east [Km]')
        plt.ylabel('north [Km]')
        ax.set_aspect('auto')
        tomo_plt = np.squeeze(tomo_cube[:, :, north_slice_index])
        if pixel_extent:
            tomo_plt = tomo_plt[:, pixel_extent[0] : pixel_extent[1]]
        ax = plt.subplot(3, 1, 2)
        plt.imshow(
            fun_plot(tomo_plt),
            interpolation='none',
            origin='lower',
            extent=[east_0, east_end, z_ax[0], z_ax[-1]],
            cmap='jet',
        )
        plt.title('north = {:3.0f} [Km]'.format(north_ax[north_slice_index] / 1000))
        plt.xlabel('east [Km]')
        plt.ylabel('z [m]')
        ax.set_aspect('auto')
        tomo_plt = np.squeeze(tomo_cube[:, east_slice_index, :])
        if pixel_extent:
            tomo_plt = tomo_plt[:, pixel_extent[2] : pixel_extent[3]]
        ax = plt.subplot(3, 1, 3)
        plt.imshow(
            fun_plot(tomo_plt),
            interpolation='none',
            origin='lower',
            extent=[north_0, north_end, z_ax[0], z_ax[-1]],
            cmap='jet',
        )
        plt.title('east = {:3.0f} [Km]'.format(east_ax[east_slice_index] / 1000))
        plt.xlabel('north [Km]')
        plt.ylabel('z [m]')
        ax.set_aspect('auto')
        plt.subplots_adjust(hspace=0.4, wspace=0.2)

        plt.figure()
        tomo_plt = np.squeeze(tomo_cube[:, east_slice_index, north_slice_index])
        plt.plot(fun_plot(tomo_plt))
        plt.title(
            'north = {} [Km], east = {} [Km]'.format(
                north_ax[north_slice_index] / 1000, east_ax[east_slice_index] / 1000
            )
        )
        plt.xlabel('z [m]')
        if plot_flag == 'db':
            plt.ylabel(plot_flag)
        ax.set_aspect('auto')

        if save_power_cube_path:
            np.save(os.path.join(save_power_cube_path, 'tomo_cube.npy'), tomo_cube)
