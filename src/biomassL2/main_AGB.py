# import libraries
import os
import numpy as np
import logging
import shutil
import scipy as sp
import scipy.optimize
import scipy.interpolate
import scipy.ndimage
import scipy.misc
import scipy.signal
from gdalconst import GA_ReadOnly
from osgeo import gdal
from equi7grid.equi7grid import Equi7Grid
from equi7grid.image2equi7grid import image2equi7grid
from scipy.signal import convolve2d

# biomassL2 processor imports
from biomassL2.processing_AGB import (
    tableLookupInt,
    compute_processing_blocs_order,
    forward,
    inverse,
    regularizeIndices,
    tableLookupFloat,
    interp2d_wrapper,
    check_intersection,
    merge_agb_intermediate,
    mean_on_rois,
    initialize_inversion_parameters,
    get_projection_from_path,
)
from biomassL2.utility_orchestrator import (
    read_and_oversample_data,
    read_and_oversample_aux_data,
    choose_equi7_sampling,
    check_if_path_exists,
    geocoding_init,
    geocoding,
    tiff_formatter,
    fnf_equi7_load_filter_equi7format,
    fnf_tandemx_load_filter_equi7format,
    apply_calibration_screens,
    apply_dem_flattening,
    compute_and_oversample_geometry_auxiliaries,
    check_if_geometry_auxiliaries_are_present,
    check_fnf_folder_format,
    check_cal_format,
    get_min_time_stamp_repository,
    get_raster_cal_names,
    get_foss_cal_names,
    resolution_heading_correction,
)
from biomassL2.IO_interfaces import (
    parse_chains_input_file,
    parse_chains_configuration_file,
    decode_unique_acquisition_id_string,
)
from biomassL2.utility_functions import save_breakpoints
from biomassL2.ground_cancellation import ground_cancellation


def main_AGB(input_file_xml, configuration_file_xml, geographic_boundaries, geographic_boundaries_per_stack, gdal_path):

    ########################## INITIAL STEPS ##############################

    logging.info('AGB: Reading chains configuration files')
    check_if_path_exists(configuration_file_xml, 'FILE')
    proc_conf = parse_chains_configuration_file(configuration_file_xml)
    proc_inputs = parse_chains_input_file(input_file_xml)

    ### managing output folders:
    products_folder = os.path.join(proc_inputs.output_folder, 'Products')
    if proc_conf.save_breakpoints:
        breakpoints_output_folder = os.path.join(products_folder, 'breakpoints')
        logging.info('AGB: Breakpoints will be saved into: ' + breakpoints_output_folder)
        os.makedirs(breakpoints_output_folder)

    temp_output_folder = os.path.join(products_folder, 'temporary_processing')
    logging.info('AGB: Temporary data folder:' + temp_output_folder + '\n')

    os.makedirs(temp_output_folder)

    ### get temporal date time of the input data (get the minimum date from all the stacks)
    time_tag_mjd_initial = get_min_time_stamp_repository(proc_inputs.L1c_repository, proc_inputs.stack_composition)

    # read the cals:
    cal_format = check_cal_format(proc_inputs.reference_agb_folder)
    if cal_format == 'GeoJSON':
        cal_fnames = get_foss_cal_names(proc_inputs.reference_agb)
        flag_cal = 0
        pass
    elif cal_format == 'RASTER':
        flag_cal = 1
        cal_fnames = get_raster_cal_names(proc_inputs.reference_agb_folder)

        # LUT: CAL
        lut_cal = np.zeros((len(cal_fnames), 5))
        lut_cal_paths = []
        for idx_cal, equi7_cal_fname in enumerate(cal_fnames):

            data_driver = gdal.Open(equi7_cal_fname, GA_ReadOnly)
            # [ upper_left_easting_coord, sampling_step_east_west, 0, upper_left_northing_coord, 0, sampling_step_north_south ]
            geotransform_equi7 = data_driver.GetGeoTransform()
            east_len = data_driver.RasterXSize
            north_len = data_driver.RasterYSize

            data_driver = None

            east_min = geotransform_equi7[0]
            east_max = east_min + geotransform_equi7[1] * east_len
            north_in = geotransform_equi7[3]
            north_out = north_in + geotransform_equi7[5] * north_len

            lut_cal[idx_cal, :] = [east_min, east_max, min(north_in, north_out), max(north_in, north_out), flag_cal]
            lut_cal_paths.append(equi7_cal_fname)

            logging.info('using calibration file: ' + equi7_cal_fname)

    ### initialize the equi7 sampling grid
    equi7_sampling_intermediate = choose_equi7_sampling(
        proc_conf.AGB.intermediate_ground_averaging, proc_conf.AGB.intermediate_ground_averaging / 2
    )
    e7g_intermediate = Equi7Grid(equi7_sampling_intermediate)
    logging.info('EQUI7 Grid sampling used for intermediate products: {}'.format(equi7_sampling_intermediate))

    # it is loaded if equi7, or converted to equi7 if TANDEM-X
    # it is a list containing all the loaded FNF-FTILES
    fnf_format = check_fnf_folder_format(proc_inputs.forest_mask_catalogue_folder)
    if fnf_format == 'TANDEM-X':
        logging.info('Initial Forest mask is in TANDEM-X format, converting to equi7...')

        # conversion step 1: get the tiff names of current zone
        equi7_fnf_mask_fnames = fnf_tandemx_load_filter_equi7format(
            proc_inputs.forest_mask_catalogue_folder,
            e7g_intermediate,
            proc_conf.AGB.intermediate_ground_averaging,
            temp_output_folder,
            gdal_path,
            geographic_boundaries,
            time_tag_mjd_initial,
        )

    elif fnf_format == 'EQUI7':
        logging.info('Initial Forest mask reading and formatting...')

        # re format the input mask in equi7 with the output resolution:
        equi7_fnf_mask_fnames = fnf_equi7_load_filter_equi7format(
            proc_inputs.forest_mask_catalogue_folder,
            e7g_intermediate,
            proc_conf.AGB.intermediate_ground_averaging,
            temp_output_folder,
            gdal_path,
        )

    # LUT: FNF
    lut_fnf = np.zeros((len(equi7_fnf_mask_fnames), 4))
    lut_fnf_paths = []
    for fnf_idx, fnf_name_curr in enumerate(equi7_fnf_mask_fnames):

        data_driver = gdal.Open(fnf_name_curr, GA_ReadOnly)
        geotransform_equi7 = data_driver.GetGeoTransform()
        east_len = data_driver.RasterXSize
        north_len = data_driver.RasterYSize

        data_driver = None

        east_min = geotransform_equi7[0]
        east_max = east_min + geotransform_equi7[1] * east_len
        north_in = geotransform_equi7[3]
        north_out = north_in + geotransform_equi7[5] * north_len

        lut_fnf[fnf_idx, :] = [east_min, east_max, min(north_in, north_out), max(north_in, north_out)]
        lut_fnf_paths.append(fnf_name_curr)

    ########################## INITIAL STEPS END #############################

    ########################## STACK BASED STEPS ##############################

    # LookUp tables initialization
    number_of_stacks = len(proc_inputs.stack_composition.keys())
    lut_progressive_stacks = np.zeros((number_of_stacks, 6))
    lut_stacks_boundaries = np.zeros((number_of_stacks, 4))
    lut_stacks_paths = []

    # sigma0 and theta initialization
    sigma0_equi7_fnames = {}
    theta_equi7_fnames = {}
    # cycle for each stack
    # for stack_key, scene_dict in proc_inputs.stacks_scenes_fnames_orders.items():
    for stack_progressive_idx, (unique_stack_id, acquisitions_pf_names) in enumerate(
        proc_inputs.stack_composition.items()
    ):

        (
            global_cycle_idx,
            heading_deg,
            rg_swath_idx,
            rg_sub_swath_idx,
            az_swath_idx,
            _,
        ) = decode_unique_acquisition_id_string(unique_stack_id + '_BSL_00')

        # LUT: progresive stacks
        lut_progressive_stacks[stack_progressive_idx, :] = [
            stack_progressive_idx,
            global_cycle_idx,
            heading_deg,
            rg_swath_idx,
            rg_sub_swath_idx,
            az_swath_idx,
        ]

        # make temporary sub-directories
        temp_output_folder_gr = os.path.join(temp_output_folder, 'ground_range_geometry', unique_stack_id)
        temp_output_folder_e7 = os.path.join(temp_output_folder, 'ground_equi7_geometry', unique_stack_id)
        os.makedirs(temp_output_folder_gr)
        os.makedirs(temp_output_folder_e7)

        ### load data ( and oversample if requested and if needed )
        try:
            logging.info('AGB: Data loading for stack ' + unique_stack_id + '; this may take a while:')

            beta0_calibrated, master_id, raster_info, raster_info_orig = read_and_oversample_data(
                proc_inputs.L1c_repository, acquisitions_pf_names, proc_conf.enable_resampling
            )

        except Exception as e:
            logging.error('AGB: error during input data reading: ' + str(e), exc_info=True)
            raise

        ### load or compute auxiliary data
        try:

            read_ref_h = not proc_conf.apply_calibration_screen and proc_conf.DEM_flattening
            read_cal_screens = proc_conf.apply_calibration_screen
            geometry_aux_are_present = check_if_geometry_auxiliaries_are_present(
                proc_inputs, unique_stack_id, acquisitions_pf_names, read_ref_h=read_ref_h
            )

            if proc_conf.compute_geometry or not geometry_aux_are_present:

                # messages for the log:
                if proc_conf.compute_geometry:
                    logging.info('AGB: calling geometry library for stack ' + unique_stack_id + '\n')
                    if geometry_aux_are_present:
                        logging.warning('    geometry auxiliaries will be overwritten for stack ' + unique_stack_id)
                    else:
                        logging.info('\n')
                else:
                    logging.warning(
                        'AGB: calling geometry library since AuxiliaryProductsFolder "Geometry" is empty or not complete \n'
                    )

                _, _, ellipsoid_slope, _, _, _, sar_geometry_master = compute_and_oversample_geometry_auxiliaries(
                    proc_inputs.L1c_repository,
                    proc_inputs,
                    unique_stack_id,
                    acquisitions_pf_names,
                    master_id,
                    proc_conf.enable_resampling,
                    force_ellipsoid=True,
                )

                (
                    ecef_grid,
                    off_nadir_angle_rad,
                    slope,
                    kz,
                    reference_height,
                    _,
                    _,
                ) = compute_and_oversample_geometry_auxiliaries(
                    proc_inputs.L1c_repository,
                    proc_inputs,
                    unique_stack_id,
                    acquisitions_pf_names,
                    master_id,
                    proc_conf.enable_resampling,
                    comp_ref_h=read_ref_h,
                    sar_geometry_master=sar_geometry_master,
                )

                logging.info('Geometry library: correcting geometry global / local reference')
                slope = slope - ellipsoid_slope
                for swath_id in kz.keys():
                    kz[swath_id] = (
                        kz[swath_id]
                        * np.sin(off_nadir_angle_rad[master_id])
                        / np.sin(off_nadir_angle_rad[master_id] - ellipsoid_slope)
                    )
                off_nadir_angle_rad[master_id] = off_nadir_angle_rad[master_id] - ellipsoid_slope
                del ellipsoid_slope

                logging.info('AGB: ...geometry auxiliaries computation done.')

            else:

                logging.info(
                    'AGB: geometry auxiliaries are provided from user, so they are now being loaded and not computed, for stack '
                    + unique_stack_id
                    + '\n'
                )

                (
                    ecef_grid,
                    off_nadir_angle_rad,
                    slope,
                    kz,
                    reference_height,
                    _,
                    _,
                    _,
                    _,
                    _,
                    _,
                    _,
                ) = read_and_oversample_aux_data(
                    proc_inputs,
                    unique_stack_id,
                    acquisitions_pf_names,
                    proc_conf.enable_resampling,
                    raster_info_orig,
                    read_ref_h=read_ref_h,
                )
                logging.info('AGB: ...geometry auxiliaries loading done.')

            # read the rest of auxiliaries which are notpart of the geometry library:
            if read_cal_screens:

                logging.warning('AGB: loading calibration screens \n')

                _, _, _, _, _, _, _, cal_screens, cal_screens_raster_info, _, _, _ = read_and_oversample_aux_data(
                    proc_inputs,
                    unique_stack_id,
                    acquisitions_pf_names,
                    proc_conf.enable_resampling,
                    raster_info_orig,
                    read_cal_screens=read_cal_screens,
                    read_ecef=False,
                    read_off_nadir=False,
                    read_slope=False,
                    read_kz=False,
                    read_ref_h=False,
                    read_dist=False,
                )
                logging.info('...done')

        except Exception as e:
            logging.error('AGB: error during auxiliary data computation and/or loading: ' + str(e), exc_info=True)
            raise

        ### Screen calibration (ground steering)
        try:
            if proc_conf.apply_calibration_screen:
                logging.info('AGB: applying calibration screen...')
                beta0_calibrated = apply_calibration_screens(
                    beta0_calibrated, raster_info, cal_screens, cal_screens_raster_info, master_id
                )
                logging.info('...done.\n')

            elif proc_conf.DEM_flattening:
                logging.info('AGB: DEM flattening... ')
                beta0_calibrated = apply_dem_flattening(beta0_calibrated, kz, reference_height, master_id, raster_info)
                logging.info('...done.\n')

        except Exception as e:
            logging.error('AGB: error during screen calibration or DEM flattening.' + str(e), exc_info=True)
            raise

        ### ground notching
        try:

            logging.info('AGB: ground contribute cancellation...:')

            DN_beta0_notched = ground_cancellation(
                beta0_calibrated,
                kz,
                proc_conf.ground_cancellation.multi_master_flag,
                proc_conf.ground_cancellation.enhanced_forest_height,
                proc_conf.ground_cancellation.equalization_flag,
                raster_info.resolution_m_slant_rg,
                off_nadir_angle_rad[master_id],
                slope,
            )  # beta0_calibrated: 3 pol (nrg x Naz  x Nimmagini); DN_beta0_notched: 3 pol (Nrg x Naz x 1 immagine)

            del beta0_calibrated, kz

        except Exception as e:
            logging.error('AGB: error during ground cancellation: ' + str(e), exc_info=True)
            raise
        logging.info('...done.\n')

        ### compute mean look angle
        look_angle_rad = np.nanmean(off_nadir_angle_rad[master_id])
        logging.info('AGB: look angle used is {} [deg] \n'.format(np.rad2deg(look_angle_rad)))

        if proc_conf.AGB.intermediate_ground_averaging > proc_conf.AGB.product_resolution / 2:
            sigma_ground_res_m = proc_conf.AGB.product_resolution / 2
            logging.warning(
                '"intermediate_ground_averaging" cannot be greater than "product_resolution/2", setting it to {}'.format(
                    sigma_ground_res_m
                )
            )
        else:
            sigma_ground_res_m = proc_conf.AGB.intermediate_ground_averaging

        if proc_conf.multilook_heading_correction:
            sigma_ground_res_m = resolution_heading_correction(sigma_ground_res_m, heading_deg)

        windtm_x = np.int(np.round(sigma_ground_res_m / raster_info.pixel_spacing_az / 2) * 2 + 1)
        windtm_y = np.int(
            np.round(sigma_ground_res_m / (raster_info.pixel_spacing_slant_rg / np.sin(look_angle_rad)) / 2) * 2 + 1
        )

        sub_factor_x = np.int((windtm_x - 1) / 2)
        sub_factor_y = np.int((windtm_y - 1) / 2)

        logging.info('AGB: multilooking of incidence angle...')
        theta_multi_looked_sr = convolve2d(
            off_nadir_angle_rad[master_id] - slope, np.ones((windtm_y, windtm_x)) / windtm_y / windtm_x, mode='same'
        )
        theta_multi_looked_sr = theta_multi_looked_sr[::sub_factor_y, ::sub_factor_x]
        logging.info('...done.\n')

        sigma0_sr = {}
        for pol_name in DN_beta0_notched.keys():

            logging.info('AGB: multilooking of ground notched for polarization {}...'.format(pol_name))
            beta0_notched_multi_looked = convolve2d(
                np.absolute(DN_beta0_notched[pol_name]) ** 2,
                np.ones((windtm_y, windtm_x)) / windtm_y / windtm_x,
                mode='same',
            )

            logging.info('AGB: sigma0 computation for polarization {}...'.format(pol_name))
            sigma0_sr[pol_name] = beta0_notched_multi_looked[::sub_factor_y, ::sub_factor_x] * np.sin(
                theta_multi_looked_sr
            )

            sigma0_sr[pol_name][sigma0_sr[pol_name] < 0] = np.NaN

            logging.info('...done.\n')

        del beta0_notched_multi_looked

        ### saving breakpoints
        if proc_conf.save_breakpoints:
            logging.info('AGB: saving breakpoints (in slant range geometry) on ' + breakpoints_output_folder)
            post_string = '_SR_' + unique_stack_id

            breakpoint_names = ['ground_cancelled_data' + post_string]

            save_breakpoints(breakpoints_output_folder, breakpoint_names, [DN_beta0_notched])
            logging.info('...done.\n')

        az_vec_subs = np.arange(0, raster_info.num_lines, sub_factor_x)
        rg_vec_subs = np.arange(0, raster_info.num_samples, sub_factor_y)

        ### Interpolate it over a regular lat lon grid (with grid data): generate a regular grid for the interpolation by using Max and Min lat lon from the ECEFGRID_LLH (make the vectors a bit longer), for the grid steps use the minimum steps from the ECEFGRID_LLH)
        logging.info(unique_stack_id + ': Geocoding data...')
        try:

            # initialize the geocoding
            min_spacing_m = min(
                sub_factor_x * raster_info.pixel_spacing_az, sub_factor_y * raster_info.pixel_spacing_slant_rg
            )
            min_spacing_m = min(min_spacing_m, equi7_sampling_intermediate)

            logging.info('Geocoding spacing set to {} [m]'.format(min_spacing_m))

            (
                lon_in,
                lat_in,
                alt_in,
                lonMeshed_out,
                latMeshed_out,
                valid_values_mask,
                lon_regular_vector,
                used_lon_step,
                lat_regular_vector,
                used_lat_step,
            ) = geocoding_init(ecef_grid, rg_vec_subs, az_vec_subs, min_spacing_m)

            # geocode the sigma0 (three polarizations)
            sigma0_gr = {}
            for pol_name in sigma0_sr.keys():

                logging.info('AGB: geocoding the sigma0 for polarization {}...'.format(pol_name))
                sigma0_gr[pol_name] = geocoding(
                    sigma0_sr[pol_name], lon_in, lat_in, lonMeshed_out, latMeshed_out, valid_values_mask
                )

                logging.info('...done.\n')

            del sigma0_sr

            # geocode the theta incidence angle
            logging.info('AGB: Geocoding of incidence angle...')
            theta_multi_looked_gr = geocoding(
                theta_multi_looked_sr, lon_in, lat_in, lonMeshed_out, latMeshed_out, valid_values_mask
            )
            logging.info('...done.\n')
            del theta_multi_looked_sr

            logging.info('...done.\n')

        except Exception as e:
            logging.error('FH: error during geocoding: ' + str(e), exc_info=True)
            raise

        ### create GEOTIFF of all the layers (estimation, mask ):
        logging.info(unique_stack_id + ': formatting data to GEOTIFF...')
        try:

            upper_left_easting_coord = lon_regular_vector[0]  # i.e. horizontal
            sampling_step_east_west = used_lon_step
            upper_left_northing_coord = lat_regular_vector[0]  # i.e. vertical
            sampling_step_north_south = used_lat_step
            geotransform = [
                upper_left_easting_coord,
                sampling_step_east_west,
                0,
                upper_left_northing_coord,
                0,
                sampling_step_north_south,
            ]

            # geotiff of the sigma0 (three polarizations)
            sigma0_ground_fnames = {}
            for pol_name in sigma0_gr.keys():
                sigma0_ground_fnames[pol_name] = os.path.join(
                    temp_output_folder_gr, 'sigma0_' + pol_name + '_' + unique_stack_id + '.tif'
                )

                tiff_formatter(
                    sigma0_gr[pol_name], sigma0_ground_fnames[pol_name], geotransform, gdal_data_format=gdal.GDT_Float32
                )

            del sigma0_gr

            # geotiff of the theta
            theta_ground_fname = os.path.join(temp_output_folder_gr, 'theta' + '_' + unique_stack_id + '.tif')

            tiff_formatter(theta_multi_looked_gr, theta_ground_fname, geotransform, gdal_data_format=gdal.GDT_Float32)

            del theta_multi_looked_gr
            logging.info('...done.\n')

        except Exception as e:
            logging.error('FH: error during GEOTIFF formatting: ' + str(e), exc_info=True)
            raise

        ### formatting data to EQUI7
        logging.info(unique_stack_id + ': formatting into EQUI7 grid...')
        try:

            sigma0_equi7_fnames[unique_stack_id] = {}
            # equi7 of the sigma0 (three polarizations)
            for pol_name in sigma0_ground_fnames.keys():

                equi7_sigma0_outdir = os.path.join(temp_output_folder_e7, 'sigma0_' + pol_name)

                logging.info('image2equi7grid IN: ' + sigma0_ground_fnames[pol_name] + ' , OUT:' + equi7_sigma0_outdir)
                sigma0_equi7_fnames[unique_stack_id][pol_name] = image2equi7grid(
                    e7g_intermediate,
                    sigma0_ground_fnames[pol_name],
                    equi7_sigma0_outdir,
                    gdal_path=gdal_path,
                    inband=None,
                    subgrid_ids=None,
                    accurate_boundary=False,
                    resampling_type='bilinear',
                    tile_nodata=np.nan,
                )

            # equi7 of the theta
            equi7_theta_outdir = os.path.join(temp_output_folder_e7, 'theta')

            logging.info('image2equi7grid IN: ' + theta_ground_fname + ' , OUT:' + equi7_theta_outdir)
            theta_equi7_fnames[unique_stack_id] = image2equi7grid(
                e7g_intermediate,
                theta_ground_fname,
                equi7_theta_outdir,
                gdal_path=gdal_path,
                inband=None,
                subgrid_ids=None,
                accurate_boundary=False,
                resampling_type='bilinear',
                tile_nodata=np.nan,
            )

        except Exception as e:
            logging.error('AGB: error during EQUI7 formatting: ' + str(e), exc_info=True)
            raise

        if not os.path.exists(theta_equi7_fnames[unique_stack_id][0]):
            error_message = 'EQUI7 grid has not been generated, output is absent '
            logging.error(error_message)
            raise

        # N0: upper left corner for output map (UTM 32S)
        _, x_upper_left, y_upper_left = e7g_intermediate.lonlat2xy(
            geographic_boundaries_per_stack[unique_stack_id].lon_min,
            geographic_boundaries_per_stack[unique_stack_id].lat_max,
        )
        north_in = y_upper_left
        east_min = x_upper_left

        # lower right corner for output map (UTM 32S)
        _, x_lower_right, y_lower_right = e7g_intermediate.lonlat2xy(
            geographic_boundaries_per_stack[unique_stack_id].lon_max,
            geographic_boundaries_per_stack[unique_stack_id].lat_min,
        )
        north_out = y_lower_right
        east_max = x_lower_right

        lut_stacks_boundaries[stack_progressive_idx, :] = [
            east_min,
            east_max,
            min(north_in, north_out),
            max(north_in, north_out),
        ]

        lut_stacks_paths.append(os.path.join(temp_output_folder_e7))

        logging.info('...done.\n')

    ########################## END OF STACK BASED STEPS ######################

    ########################## AGB INVERSION CORE ############################

    ### INITIALIZE AGB INVERSION ###

    # output folders:
    output_folder = os.path.join(products_folder, 'global_AGB')
    os.makedirs(output_folder)
    temp_output_folder_agb_est = os.path.join(temp_output_folder, 'agb_estimation', unique_stack_id)
    os.makedirs(temp_output_folder_agb_est)

    sigma_tab_pixel_folder = os.path.join(temp_output_folder_agb_est, 'sigma_tab')
    os.makedirs(sigma_tab_pixel_folder)
    theta_tab_pixel_folder = os.path.join(temp_output_folder_agb_est, 'theta_tab')
    os.makedirs(theta_tab_pixel_folder)
    acq_tab_pixel_folder = os.path.join(temp_output_folder_agb_est, 'acq_tab')
    os.makedirs(acq_tab_pixel_folder)
    polid_tab_pixel_folder = os.path.join(temp_output_folder_agb_est, 'polid_tab')
    os.makedirs(polid_tab_pixel_folder)

    # parameter text file path
    par_path = os.path.join(output_folder, 'par_est.txt')

    # read and initialize all the parameters needed for the inversion
    (
        dx_in,
        N0,
        E0,
        NN,
        EN,
        dE,
        dN,
        dE_roi,
        dN_roi,
        wE_roi,
        wN_roi,
        dE_par,
        dN_par,
        wE_par,
        wN_par,
        all_vars,
        l_lims,
        a_lims,
        n_lims,
        w_lims,
        N_tests,
        geographic_grid_sampling,
        sub_grid_string,
    ) = initialize_inversion_parameters(
        equi7_sampling_intermediate, proc_inputs.geographic_grid_sampling, geographic_boundaries, proc_conf.AGB
    )

    ### initialize the equi7 sampling grid or output product
    e7g_product = Equi7Grid(geographic_grid_sampling)
    logging.info('EQUI7 Grid sampling used for final products: {}'.format(geographic_grid_sampling))

    equi7_subgrid_names_list = []
    for lut_stacks_path in lut_stacks_paths:
        for equi7_subgrid_name_curr in os.listdir(os.path.join(lut_stacks_path, 'sigma0_hh')):

            equi7_subgrid_names_list.append(equi7_subgrid_name_curr)

    if len(set(equi7_subgrid_names_list)) > 1:
        err_str = 'AGB: biomass L2 prototype does not support multiple equi7 sub grid regions'
        logging.error(err_str)
        raise
    subgrid_code = equi7_subgrid_names_list[0][6:8]
    projection_prev = ''
    equi7_agb_est_out_tif_names_not_merged = []

    ### INITIALIZE AGB INVERSION END ###

    ### PREPARE INVERSION MODULE ###

    # parameter block axis (monidimensional axis, coordinates and bidimensional meshes)
    EE_par_mesh, NN_par_mesh = np.meshgrid(np.arange(E0, EN, dE_par), np.arange(N0, NN, dN_par))
    EE_par_coordinates, NN_par_coordinates = [np.float64(x.flatten()) for x in [EE_par_mesh, NN_par_mesh]]
    # if just one block is present, recast the coordinates as numopy array of one element
    if not EE_par_coordinates.shape:
        EE_par_coordinates = np.array([EE_par_coordinates])
        NN_par_coordinates = np.array([NN_par_coordinates])

    # derived quantities
    dl, da, dn, dw = [np.diff(lims)[0] for lims in [l_lims, a_lims, n_lims, w_lims]]
    l0, a0, n0, w0 = [lims[0] for lims in [l_lims, a_lims, n_lims, w_lims]]

    # read acquisition info table
    acq_info_tab = lut_progressive_stacks
    # restructure so that polarisation also is included (replicate rows three times)
    acq_info_tab = np.column_stack(
        (np.kron(np.arange(3), np.ones(acq_info_tab.shape[0])), np.kron(np.ones((3, 1)), acq_info_tab))
    )

    # Compute the starting block and the nblock processing order
    par_block_index_curr, block_order = compute_processing_blocs_order(
        lut_cal, EE_par_coordinates, wE_par, NN_par_coordinates, wN_par
    )

    ### PREPARE INVERSION MODULE END ###

    ### STARTING THE MAIN INVERSION ALGORITHM ###

    # set up "finished" table for all parameter blocks
    number_of_steps = len(EE_par_coordinates.flatten())
    par_finished = np.zeros((number_of_steps), dtype='bool')
    # run as long as the current block has a positive value
    counter_iteration = 0
    while ~np.all(par_finished):

        # show status
        counter_iteration = counter_iteration + 1
        logging.info('AGB: Running block #{}'.format(par_block_index_curr))
        logging.info('AGB: estimation step {} of {}'.format(counter_iteration, number_of_steps))

        ## SELECT VALID PIXELS AND ROIS FOR CURRENT BLOCK
        # extent of current block
        par_extent_curr = np.array(
            [
                EE_par_coordinates[par_block_index_curr],
                EE_par_coordinates[par_block_index_curr] + wE_par,
                NN_par_coordinates[par_block_index_curr],
                NN_par_coordinates[par_block_index_curr] + wN_par,
            ]
        )  # [min_east, max_east, min_north, max_north]

        # Generate roi grid for current parameter block
        EE_roi_mesh, NN_roi_mesh = np.meshgrid(
            np.arange(par_extent_curr[0], par_extent_curr[1], dE_roi),
            np.arange(par_extent_curr[2], par_extent_curr[3], dN_roi),
        )

        # ROI axis for current parameter block
        EE_roi_axis = np.arange(par_extent_curr[0], par_extent_curr[1], dE_roi)
        NN_roi_axis = np.arange(par_extent_curr[2], par_extent_curr[3], dN_roi)

        # Pixel axis for current parameter block)
        EE_pixel_axis = np.arange(par_extent_curr[0], par_extent_curr[1], dE)
        NN_pixel_axis = np.arange(par_extent_curr[2], par_extent_curr[3], dN)

        EE_pixel_mesh, NN_pixel_mesh = np.meshgrid(EE_pixel_axis, NN_pixel_axis)

        number_of_rois = len(EE_roi_axis) * len(NN_roi_axis)
        if number_of_rois < proc_conf.AGB.min_number_of_rois:
            logging.info(
                'skipping block #{} because the number of rois #{} cannot be less than #{}'.format(
                    par_block_index_curr, number_of_rois, proc_conf.AGB.min_number_of_rois
                )
            )
            # swap the flag
            par_finished[par_block_index_curr] = True
            # remove current par block from list
            block_order = block_order[block_order != par_block_index_curr]
            # select next par block
            if len(block_order) > 0:
                # if there is at least one left, just take the next closest to CALdata
                par_block_index_curr = block_order[0]
                continue
            else:
                break

        EE_pixel_mesh, NN_pixel_mesh = np.meshgrid(EE_pixel_axis, NN_pixel_axis)

        try:

            # Initialize tables creation cycles
            number_of_stacks = len(lut_stacks_paths)
            block_has_data = np.zeros(number_of_stacks, dtype='bool')
            for stack_idx, stack_path in enumerate(lut_stacks_paths):

                # go ahead only if current parameter block is at least partially contained in the data stack:
                block_has_data[stack_idx] = check_intersection(
                    lut_stacks_boundaries[stack_idx, 0],
                    lut_stacks_boundaries[stack_idx, 1],
                    lut_stacks_boundaries[stack_idx, 2],
                    lut_stacks_boundaries[stack_idx, 3],
                    par_extent_curr[0],
                    par_extent_curr[1],
                    par_extent_curr[3],
                    par_extent_curr[2],
                )

            block_has_cal = np.zeros(len(lut_cal_paths), dtype='bool')
            for cal_idx, cal_path in enumerate(lut_cal_paths):

                # go ahead only if current parameter block is at least partially contained in the data stack:
                block_has_cal[cal_idx] = check_intersection(
                    lut_cal[cal_idx, 0],
                    lut_cal[cal_idx, 1],
                    lut_cal[cal_idx, 2],
                    lut_cal[cal_idx, 3],
                    par_extent_curr[0],
                    par_extent_curr[1],
                    par_extent_curr[3],
                    par_extent_curr[2],
                )

            number_of_stacks_inside = np.sum(block_has_data == True)

            if not number_of_stacks_inside:
                logging.info('skipping block #{} due to too invalid data points'.format(par_block_index_curr))
                # swap the flag
                par_finished[par_block_index_curr] = True
                # remove current par block from list
                block_order = block_order[block_order != par_block_index_curr]
                # select next par block
                if len(block_order) > 0:
                    # if there is at least one left, just take the next closest to CALdata
                    par_block_index_curr = block_order[0]
                    continue
                else:
                    break

            equi7_grid_name = os.listdir(os.path.join(lut_stacks_paths[0], 'sigma0_hh'))[0]

            pol_names = ['hh', 'vh', 'vv']

            acqid_tab = np.nan * np.zeros(
                (number_of_stacks_inside * len(pol_names), number_of_rois)
            )  # progressive stacks ids (replicates values for each roi)
            sigma_tab = np.nan * np.zeros(
                (number_of_stacks_inside * len(pol_names), number_of_rois)
            )  # sigma0 values for each roi and for each pol in each stack
            theta_tab = np.nan * np.zeros(
                (number_of_stacks_inside, number_of_rois)
            )  # theta values for each roi in each stack: it is replicated after for matching the three polarizations shape

            number_of_pixels = len(EE_pixel_axis) * len(NN_pixel_axis)
            sigma_tab_pixels = np.nan * np.zeros((number_of_stacks_inside * len(pol_names), number_of_pixels))
            theta_tab_pixels = np.nan * np.zeros((number_of_stacks_inside, number_of_pixels))
            acqid_tab_pixels = np.nan * np.zeros((number_of_stacks_inside * len(pol_names), number_of_pixels))

            # loading all fnf masks that fall inside a block: note, more than ine fnf equi7 tile can fall inside
            # loaded masks needs to be re-interpolated over the pixels roi grid (EE_pixel_axis, NN_pixel_axis)
            fnf_mask_data_interp = np.zeros((len(NN_pixel_axis), len(EE_pixel_axis)), dtype='float')

            counter_fnf = 0
            for fnf_idx, fnf_path in enumerate(lut_fnf_paths):

                fnf_boundaries = lut_fnf[fnf_idx]

                fnf_is_inside = check_intersection(
                    fnf_boundaries[0],
                    fnf_boundaries[1],
                    fnf_boundaries[2],
                    fnf_boundaries[3],
                    par_extent_curr[0],
                    par_extent_curr[1],
                    par_extent_curr[3],
                    par_extent_curr[2],
                )
                if fnf_is_inside:

                    fnf_mask_data_interp_curr = np.round(
                        interp2d_wrapper(fnf_path, 1, EE_pixel_axis, NN_pixel_axis, fill_value=float(0))
                    )

                    # mean all the fnf tiles
                    fnf_mask_data_interp = np.ceil(
                        merge_agb_intermediate(fnf_mask_data_interp, fnf_mask_data_interp_curr, method='nan_mean')
                    )

                    fnf_mask_data_interp[fnf_mask_data_interp != 1] = 0
                    fnf_mask_data_interp[np.isnan(fnf_mask_data_interp)] = 0

                    counter_fnf = counter_fnf + 1

            if counter_fnf == 0:

                err_str = 'Cannot find any FNF mask falling in current block coordinates.'
                logging.error(err_str)
                raise ImportError(err_str)

            # cycle for each stack
            counter_stacks = 0
            counter_stack_pols = 0
            for stack_idx, stack_path in enumerate(lut_stacks_paths):

                # go ahead only if current stack is (at least partially) contained in the current parameter block:
                if block_has_data[stack_idx]:

                    progressive_stack_idx = int(lut_progressive_stacks[stack_idx, 0])

                    # cycle for each polarization inside the stack
                    for pol_name in pol_names:

                        sigma0_eq7_parent_dir_name = os.path.join(stack_path, 'sigma0_' + pol_name, equi7_grid_name)
                        equi7_tiles_names = os.listdir(sigma0_eq7_parent_dir_name)

                        # cycle over all the equi7 tiles, interpolate over pixel grid and mean them togheter
                        stack_data_interp = np.NaN * np.zeros((len(NN_pixel_axis), len(EE_pixel_axis)), dtype='float')
                        for equi7_tile_name in equi7_tiles_names:

                            sigma0_tiff_name = os.path.join(
                                sigma0_eq7_parent_dir_name,
                                equi7_tile_name,
                                'sigma0_'
                                + pol_name
                                + '_'
                                + os.path.basename(stack_path)
                                + '_'
                                + equi7_grid_name[6:]
                                + '_'
                                + equi7_tile_name
                                + '.tif',
                            )

                            stack_data_interp_curr = interp2d_wrapper(
                                sigma0_tiff_name, 1, EE_pixel_axis, NN_pixel_axis, fill_value=np.NaN
                            )

                            stack_data_interp = merge_agb_intermediate(
                                stack_data_interp, stack_data_interp_curr, method='nan_mean'
                            )

                        # masking the stack:
                        stack_data_interp[fnf_mask_data_interp != 1] = np.NaN

                        sigma_tab_pixels[counter_stack_pols, :] = stack_data_interp.flatten()

                        # Mean on ROI
                        data_roi_means_vec = mean_on_rois(
                            stack_data_interp,
                            EE_pixel_mesh,
                            NN_pixel_mesh,
                            EE_roi_axis,
                            dE_roi,
                            NN_roi_axis,
                            dN_roi,
                            'mean',
                        )

                        # Insert in table
                        sigma_tab[counter_stack_pols, :] = data_roi_means_vec

                        progressive_stack_idx_vec = np.ones(number_of_rois, dtype='int') * progressive_stack_idx
                        progressive_stack_idx_vec[np.isnan(data_roi_means_vec)] = -1  # invalid value

                        acqid_tab[counter_stack_pols, :] = progressive_stack_idx_vec

                        progressive_stack_idx_vec = np.ones(number_of_pixels, dtype='int') * progressive_stack_idx
                        progressive_stack_idx_vec[np.isnan(stack_data_interp.flatten())] = -1  # invalid value

                        acqid_tab_pixels[counter_stack_pols, :] = progressive_stack_idx_vec

                        counter_stack_pols = counter_stack_pols + 1

                    # cycle over all the equi7 tiles, interpolate over pixel grid and mean them togheter
                    theta_eq7_parent_dir_name = os.path.join(stack_path, 'theta', equi7_grid_name)

                    stack_theta_interp = np.NaN * np.zeros((len(NN_pixel_axis), len(EE_pixel_axis)), dtype='float')
                    for equi7_tile_name in equi7_tiles_names:

                        theta_tiff_name = os.path.join(
                            theta_eq7_parent_dir_name,
                            equi7_tile_name,
                            'theta_'
                            + os.path.basename(stack_path)
                            + '_'
                            + equi7_grid_name[6:]
                            + '_'
                            + equi7_tile_name
                            + '.tif',
                        )

                        stack_theta_interp_curr = interp2d_wrapper(
                            theta_tiff_name, 1, EE_pixel_axis, NN_pixel_axis, fill_value=np.NaN
                        )

                        # mean all the equi7 tiles
                        stack_theta_interp = merge_agb_intermediate(
                            stack_theta_interp, stack_theta_interp_curr, method='nan_mean'
                        )

                        projection_curr = get_projection_from_path(theta_tiff_name)
                        if projection_prev and (projection_prev != projection_curr):
                            err_str = 'multiple equi7 projections are not supported.'
                            logging.error(err_str)
                            raise ImportError(err_str)
                        projection_prev = projection_curr

                    # masking the stack:
                    stack_theta_interp[fnf_mask_data_interp != 1] = np.NaN

                    theta_tab_pixels[counter_stacks, :] = stack_theta_interp.flatten()

                    # Mean on ROI
                    theta_roi_means_vec = mean_on_rois(
                        stack_theta_interp,
                        EE_pixel_mesh,
                        NN_pixel_mesh,
                        EE_roi_axis,
                        dE_roi,
                        NN_roi_axis,
                        dN_roi,
                        'mean',
                    )

                    # Insert in table
                    theta_tab[counter_stacks, :] = theta_roi_means_vec

                    counter_stacks = counter_stacks + 1

            # replicate equal values for all the 3 polarizations:
            theta_tab = np.kron(theta_tab, np.ones((3, 1)))
            theta_tab_pixels = np.kron(theta_tab_pixels, np.ones((3, 1)))

            N_laypol = sigma_tab.shape[0]  # number of stacks x number of polarizations'
            N_lay = np.int32(N_laypol / 3)  # number of stacks
            # create polarisation ID table matching sigma data table in shape and size
            polid_tab = np.kron(
                np.ones((N_lay, 1)), np.kron(np.array([np.arange(3)]).transpose(), np.ones((1, number_of_rois)))
            )
            polid_tab_pixels = np.kron(
                np.ones((N_lay, 1)), np.kron(np.array([np.arange(3)]).transpose(), np.ones((1, number_of_pixels)))
            )

            # save and delete all the pixel tabs:
            np.save(
                os.path.join(sigma_tab_pixel_folder, 'sigma_tab_pixels_block_{}'.format(par_block_index_curr)),
                sigma_tab_pixels,
            )
            np.save(
                os.path.join(acq_tab_pixel_folder, 'acq_tab_pixels_block_{}'.format(par_block_index_curr)),
                acqid_tab_pixels,
            )
            np.save(
                os.path.join(theta_tab_pixel_folder, 'theta_tab_pixels_block_{}'.format(par_block_index_curr)),
                theta_tab_pixels,
            )
            np.save(
                os.path.join(polid_tab_pixel_folder, 'polid_tab_pixels_block_{}'.format(par_block_index_curr)),
                polid_tab_pixels,
            )
            del sigma_tab_pixels, acqid_tab_pixels, theta_tab_pixels, polid_tab_pixels

            ## FILTER OUT INVALID ACQUISITIONS
            # find independent measurements within current tile and count them
            unique_acqs, acq_count = np.unique(acqid_tab, return_counts=True)
            # extract acquisitions meeting the requirements on minimal number of independent measurements
            valid_acqs = np.isin(
                acq_info_tab[:, 1], unique_acqs[acq_count >= proc_conf.AGB.min_number_of_rois_per_stack * 3]
            )
            # mask out invalid acquisitions (either too few measurements or there are nans in sigma or theta images)
            acqid_tab[
                (~np.isin(acqid_tab, acq_info_tab[valid_acqs, 1]))
                | np.isnan(sigma_tab)
                | (sigma_tab < 0)
                | np.isnan(theta_tab)
            ] = -1

            # mask out data in sigma0 and local tables
            sigma_tab[acqid_tab == -1] = np.nan
            theta_tab[acqid_tab == -1] = np.nan

        except Exception as e:
            logging.error('AGB: error during tables creation.' + str(e), exc_info=True)
            raise

        # make sure we have enough data
        if not sum(valid_acqs == True):
            logging.info('skipping block #{} due to too invalid data points'.format(par_block_index_curr))
            # swap the flag
            par_finished[par_block_index_curr] = True
            # remove current par block from list
            block_order = block_order[block_order != par_block_index_curr]
            # select next par block
            if len(block_order) > 0:
                # if there is at least one left, just take the next closest to CALdata
                par_block_index_curr = block_order[0]
                continue
            else:
                break

        ## PREPARE PARAMETER ID TABLES
        # create roi id table
        roiid_tab = np.int32(np.array([np.arange(number_of_rois)]) * np.ones((N_laypol, 1)))
        # create look-up tables for acquisition ID and unique parameter ID
        par_id_luts = [
            np.unique(acq_info_tab[valid_acqs, :] * np.array([curr_var]), axis=0, return_inverse=True)[1]
            for curr_var in all_vars
        ]
        # create parameter id tables matching data tables
        # acq_info_tab is ordered differently from other tables with respect to acquisition id and polarization
        # tableLookupInt function takes that into account
        par_tabs = [
            tableLookupInt(
                acq_info_tab[valid_acqs, :2], curr_lut, np.column_stack((polid_tab.flatten(), acqid_tab.flatten()))
            ).reshape(acqid_tab.shape)
            for curr_lut in par_id_luts
        ]
        # extract the number of unique values for each parameter
        all_N_par = np.int32(np.array([np.nanmax(par_tab) + 1 for par_tab in par_tabs + [roiid_tab]]))
        # calculate the offset from first parameter in the large vector
        par_offsets = np.int32(np.concatenate((np.zeros(1), np.cumsum(all_N_par))))[:-1]

        # if at least one tile has been finished, load the low-resolution estimates for the relevant pixels
        counter_cal = 0
        cal_data_interp = np.nan * np.zeros((len(NN_pixel_axis), len(EE_pixel_axis)), dtype='float')
        cal_std1_interp = np.nan * np.zeros((len(NN_pixel_axis), len(EE_pixel_axis)), dtype='float')
        cal_std2_interp = np.nan * np.zeros((len(NN_pixel_axis), len(EE_pixel_axis)), dtype='float')
        cal_std3_interp = np.nan * np.zeros((len(NN_pixel_axis), len(EE_pixel_axis)), dtype='float')
        for cal_idx, cal_path in enumerate(lut_cal_paths):
            logging.info('current block uses cal: ' + cal_path)
            if block_has_cal[cal_idx]:
                logging.info('above cal is in the block')
                counter_cal = counter_cal + 1

                cal_data_interp_curr = interp2d_wrapper(cal_path, 1, EE_pixel_axis, NN_pixel_axis, fill_value=np.NaN)
                cal_std1_interp_curr = interp2d_wrapper(cal_path, 2, EE_pixel_axis, NN_pixel_axis, fill_value=np.NaN)
                cal_std2_interp_curr = interp2d_wrapper(cal_path, 3, EE_pixel_axis, NN_pixel_axis, fill_value=np.NaN)
                cal_std3_interp_curr = interp2d_wrapper(cal_path, 4, EE_pixel_axis, NN_pixel_axis, fill_value=np.NaN)

                # mean all the equi7 tiles
                cal_data_interp = merge_agb_intermediate(cal_data_interp, cal_data_interp_curr, method='nan_mean')
                cal_std1_interp = np.sqrt(
                    merge_agb_intermediate(cal_std1_interp ** 2, cal_std1_interp_curr ** 2, method='nan_mean')
                )
                cal_std2_interp = np.sqrt(
                    merge_agb_intermediate(cal_std2_interp ** 2, cal_std2_interp_curr ** 2, method='nan_mean')
                )
                cal_std3_interp = np.sqrt(
                    merge_agb_intermediate(cal_std3_interp ** 2, cal_std3_interp_curr ** 2, method='nan_mean')
                )

        # Mean on ROI
        cal_data_roi_means_vec = mean_on_rois(
            cal_data_interp, EE_pixel_mesh, NN_pixel_mesh, EE_roi_axis, dE_roi, NN_roi_axis, dN_roi, method='nan_mean'
        )
        cal_std1_roi_means_vec = np.sqrt(
            mean_on_rois(
                cal_std1_interp ** 2,
                EE_pixel_mesh,
                NN_pixel_mesh,
                EE_roi_axis,
                dE_roi,
                NN_roi_axis,
                dN_roi,
                method='nan_mean',
            )
        )
        cal_std2_roi_means_vec = np.sqrt(
            mean_on_rois(
                cal_std2_interp ** 2,
                EE_pixel_mesh,
                NN_pixel_mesh,
                EE_roi_axis,
                dE_roi,
                NN_roi_axis,
                dN_roi,
                method='nan_mean',
            )
        )
        cal_std3_roi_means_vec = np.sqrt(
            mean_on_rois(
                cal_std3_interp ** 2,
                EE_pixel_mesh,
                NN_pixel_mesh,
                EE_roi_axis,
                dE_roi,
                NN_roi_axis,
                dN_roi,
                method='nan_mean',
            )
        )

        # replicate N_laypol times
        cal_tab = np.kron(cal_data_roi_means_vec, np.ones((N_laypol, 1)))
        std1_tab_db = np.kron(cal_std1_roi_means_vec, np.ones((N_laypol, 1)))
        std2_tab_db = np.kron(cal_std2_roi_means_vec, np.ones((N_laypol, 1)))
        std3_tab_db = np.kron(cal_std3_roi_means_vec, np.ones((N_laypol, 1)))

        ## PREPARE CALIBRATION SETS
        # now, lets take out calibration ROIs and randomise two of them to be used for calibration in multiple tests
        all_cal_ids = np.unique(roiid_tab[~np.isnan(cal_tab) & (acqid_tab > -1)])
        all_roi_ids = np.setdiff1d(np.unique(roiid_tab[acqid_tab > -1]), all_cal_ids)  # tolgo le cal dalle roi
        # number of CALs in a subset (approximately half the total number of CALs)
        N_cal_sub = np.int32(np.floor(proc_conf.AGB.fraction_of_cal_per_test / 100 * len(all_cal_ids)))
        N_roi_sub = np.int32(np.floor(proc_conf.AGB.fraction_of_roi_per_test / 100 * len(all_roi_ids)))

        logging.info('Number of cals in current block :{}'.format(N_cal_sub))
        logging.info('Number of rois in current block :{}'.format(N_roi_sub))

        # make sure we have enough data
        if not (
            (N_cal_sub >= proc_conf.AGB.min_number_of_cals_per_test)
            & (N_roi_sub >= proc_conf.AGB.min_number_of_rois_per_test)
        ):
            logging.info('skipping block #{} due to too few data points'.format(par_block_index_curr))
            # swap the flag
            par_finished[par_block_index_curr] = True
            # remove current par block from list
            block_order = block_order[block_order != par_block_index_curr]
            # select next par block
            if len(block_order) > 0:
                # if there is at least one left, just take the next closest to CALdata
                par_block_index_curr = block_order[0]
                continue
            else:
                break

        # create random subsets of N_cal and N_roi
        cal_sets = all_cal_ids[
            np.unique(
                np.row_stack(
                    [np.sort(np.random.permutation(np.arange(len(all_cal_ids)))[:N_cal_sub]) for ii in np.arange(10000)]
                ),
                axis=0,
            )
        ]
        cal_sets = cal_sets[np.random.permutation(np.arange(cal_sets.shape[0]))[:N_tests], :]
        roi_sets = all_roi_ids[
            np.unique(
                np.row_stack(
                    [np.sort(np.random.permutation(np.arange(len(all_roi_ids)))[:N_roi_sub]) for ii in np.arange(10000)]
                ),
                axis=0,
            )
        ]
        roi_sets = roi_sets[np.random.permutation(np.arange(roi_sets.shape[0]))[:N_tests], :]
        # allocate parameter table
        curr_output_table = np.nan * np.zeros((N_laypol, number_of_rois, 6, N_tests))
        # counter for tests run
        i_test = 0

        ## RUN INVERSIONS
        for cal_set, roi_set in zip(cal_sets, roi_sets):

            # create masks for valid cal and roi data
            curr_calmask = (acqid_tab != -1) & np.isin(roiid_tab, cal_set)
            curr_roimask = (acqid_tab != -1) & np.isin(roiid_tab, roi_set)
            # create zero-mean, unity-std errors of three types:
            # 1) ROI-to-ROI error
            # 2) test-to-test error
            # 3) layer-to-layer error
            eps1 = np.random.randn(1, number_of_rois)
            eps2 = np.random.randn(N_laypol, 1)
            eps3 = np.random.randn()
            # extract data for rois and cals
            s_cal = forward(sigma_tab[curr_calmask])
            c_cal = forward(np.cos(theta_tab[curr_calmask]))
            w_cal = (
                forward(np.maximum(1, cal_tab[curr_calmask]))
                + (std1_tab_db * eps1)[curr_calmask]
                + (std2_tab_db * eps2)[curr_calmask]
                + (std3_tab_db * eps3)[curr_calmask]
            )
            s_roi = forward(sigma_tab[curr_roimask])
            c_roi = forward(np.cos(theta_tab[curr_roimask]))
            # number of CAL and ROI measurements
            N = len(s_roi)
            M = len(s_cal)
            # create index vectors for rois and cals
            i0_l_cal, i0_a_cal, i0_n_cal = [
                par_tab[curr_calmask] + par_offset for par_tab, par_offset in zip(par_tabs, par_offsets)
            ]
            i0_l_roi, i0_a_roi, i0_n_roi = [
                par_tab[curr_roimask] + par_offset for par_tab, par_offset in zip(par_tabs, par_offsets)
            ]
            i0_w_roi = roiid_tab[curr_roimask] + par_offsets[-1]
            # regularize the indices to get minimal size x-vector
            i_l_cal, i_l_roi, lut_l = regularizeIndices(i0_l_cal, i0_l_roi, 0)
            i_a_cal, i_a_roi, lut_a = regularizeIndices(i0_a_cal, i0_a_roi, np.max(lut_l[:, 1]) + 1)
            i_n_cal, i_n_roi, lut_n = regularizeIndices(i0_n_cal, i0_n_roi, np.max(lut_a[:, 1]) + 1)
            i_w_roi, i_w_roi, lut_w = regularizeIndices(i0_w_roi, i0_w_roi, np.max(lut_n[:, 1]) + 1)
            # total number of parameters to be estimated
            N_par = np.max(i_w_roi) + 1

            # cost function
            def F(x):
                # calculate parameter values from x values
                l_cal = l0 + dl * np.sin(x[i_l_cal]) ** 2
                a_cal = a0 + da * np.sin(x[i_a_cal]) ** 2
                n_cal = n0 + dn * np.sin(x[i_n_cal]) ** 2
                l_roi = l0 + dl * np.sin(x[i_l_roi]) ** 2
                a_roi = a0 + da * np.sin(x[i_a_roi]) ** 2
                n_roi = n0 + dn * np.sin(x[i_n_roi]) ** 2
                w_roi = w0 + dw * np.sin(x[i_w_roi]) ** 2
                # final cost function
                A = (l_roi + a_roi * w_roi + n_roi * c_roi - s_roi) ** 2
                B = (l_cal + a_cal * w_cal + n_cal * c_cal - s_cal) ** 2
                return 1 / N * np.nansum(A) + 1 / M * np.nansum(B)

            # pre-calculate selection matrices
            X_roi = [
                np.column_stack(([np.float32(i_l_roi == ii) for ii in np.arange(N_par)])),
                np.column_stack(([np.float32(i_a_roi == ii) for ii in np.arange(N_par)])),
                np.column_stack(([np.float32(i_n_roi == ii) for ii in np.arange(N_par)])),
                np.column_stack(([np.float32(i_w_roi == ii) for ii in np.arange(N_par)])),
            ]
            X_cal = [
                np.column_stack(([np.float32(i_l_cal == ii) for ii in np.arange(N_par)])),
                np.column_stack(([np.float32(i_a_cal == ii) for ii in np.arange(N_par)])),
                np.column_stack(([np.float32(i_n_cal == ii) for ii in np.arange(N_par)])),
            ]
            # gradient of cost function
            def G(x):
                # calculate parameter values from x values
                l_cal = l0 + dl * np.sin(x[i_l_cal]) ** 2
                a_cal = a0 + da * np.sin(x[i_a_cal]) ** 2
                n_cal = n0 + dn * np.sin(x[i_n_cal]) ** 2
                l_roi = l0 + dl * np.sin(x[i_l_roi]) ** 2
                a_roi = a0 + da * np.sin(x[i_a_roi]) ** 2
                n_roi = n0 + dn * np.sin(x[i_n_roi]) ** 2
                w_roi = w0 + dw * np.sin(x[i_w_roi]) ** 2
                # derivatives of these parameters wrt x
                D_l_cal = dl * np.sin(2 * x[i_l_cal])
                D_a_cal = da * np.sin(2 * x[i_a_cal])
                D_n_cal = dn * np.sin(2 * x[i_n_cal])
                D_l_roi = dl * np.sin(2 * x[i_l_roi])
                D_a_roi = da * np.sin(2 * x[i_a_roi])
                D_n_roi = dn * np.sin(2 * x[i_n_roi])
                D_w_roi = dw * np.sin(2 * x[i_w_roi])
                # final gradient values
                A = l_roi + a_roi * w_roi + n_roi * c_roi - s_roi
                B = l_cal + a_cal * w_cal + n_cal * c_cal - s_cal
                return 2 * (
                    1
                    / N
                    * np.nansum(
                        np.array([A * D_l_roi]).transpose() * X_roi[0]
                        + np.array([A * D_a_roi * w_roi]).transpose() * X_roi[1]
                        + np.array([A * D_n_roi * c_roi]).transpose() * X_roi[2]
                        + np.array([A * D_w_roi * a_roi]).transpose() * X_roi[3],
                        axis=0,
                    )
                    + 1
                    / M
                    * np.nansum(
                        np.array([B * D_l_cal]).transpose() * X_cal[0]
                        + np.array([B * D_a_cal * w_cal]).transpose() * X_cal[1]
                        + np.array([B * D_n_cal * c_cal]).transpose() * X_cal[2],
                        axis=0,
                    )
                )

            # initial x parameter value (random vector)
            x0 = np.pi / 2 * np.random.rand(N_par)
            # perform the minimisation
            mod = sp.optimize.minimize(F, x0, jac=G)
            # show message
            logging.info('\tMinimisation for test %d/%d finished (%s)' % (i_test + 1, N_tests, mod.message))
            # extract parameter values
            x = mod.x
            # remove those that haven't changed from the initial values (meaning this parameter is not estimated)
            x[np.abs(x - x0) < 1e-4] = np.nan
            # save parameters to output matrix
            curr_output_table[curr_calmask, 0, i_test] = l0 + dl * np.sin(x[i_l_cal]) ** 2
            curr_output_table[curr_roimask, 0, i_test] = l0 + dl * np.sin(x[i_l_roi]) ** 2
            curr_output_table[curr_calmask, 1, i_test] = a0 + da * np.sin(x[i_a_cal]) ** 2
            curr_output_table[curr_roimask, 1, i_test] = a0 + da * np.sin(x[i_a_roi]) ** 2
            curr_output_table[curr_calmask, 2, i_test] = n0 + dn * np.sin(x[i_n_cal]) ** 2
            curr_output_table[curr_roimask, 2, i_test] = n0 + dn * np.sin(x[i_n_roi]) ** 2
            # save biomass values for rois
            curr_output_table[curr_roimask, 3, i_test] = w0 + dw * np.sin(x[i_w_roi]) ** 2
            # save numerator values
            curr_output_table[curr_roimask, 4, i_test] = (
                s_roi - (n0 + dn * np.sin(x[i_n_roi]) ** 2) * c_roi - (l0 + dl * np.sin(x[i_l_roi]) ** 2)
            ) * (a0 + da * np.sin(x[i_a_roi]) ** 2)
            curr_output_table[curr_calmask, 4, i_test] = (
                s_cal - (n0 + dn * np.sin(x[i_n_cal]) ** 2) * c_cal - (l0 + dl * np.sin(x[i_l_cal]) ** 2)
            ) * (a0 + da * np.sin(x[i_a_cal]) ** 2)
            # save current calibration data
            curr_output_table[curr_calmask, 5, i_test] = w_cal

            # increase test counter
            i_test += 1

        ## CONSOLIDATING RESULTS
        # calculate biomass estimates for each individual triplet and test using the equation
        agbest_tab_eq = inverse(
            np.array(
                [
                    np.nansum(curr_output_table[ilay * 3 + np.arange(3), :, 4, :], axis=0)
                    / np.nansum(curr_output_table[ilay * 3 + np.arange(3), :, 1, :] ** 2, axis=0)
                    for ilay in np.int32(np.arange(N_lay))
                ]
            )
        )
        # calculate the first residual statistics from the difference between these biomass estimates and the true biomass values for CALs
        res_1 = forward(agbest_tab_eq) - forward(np.array([[cal_tab[0, :]]]) * np.ones((N_tests, N_lay, 1))).transpose(
            [1, 2, 0]
        )
        agbstd_1 = np.nanstd(res_1)
        bias_1 = np.nanmean(res_1)
        # re-calculate biomass estimates for each individual triplet and test using the equation and bias estimates
        agbest_tab_eq = inverse(
            (
                np.array(
                    [
                        np.nansum(curr_output_table[ilay * 3 + np.arange(3), :, 4, :], axis=0)
                        / np.nansum(curr_output_table[ilay * 3 + np.arange(3), :, 1, :] ** 2, axis=0)
                        for ilay in np.int32(np.arange(N_lay))
                    ]
                )
                - bias_1
            )
        )
        # calculate biomass estimates averaged over triplets for each individual test using the equation
        agbest_tab_eqbar = inverse(
            (
                np.nansum(curr_output_table[:, :, 4, :], axis=0) / np.nansum(curr_output_table[:, :, 1, :] ** 2, axis=0)
                - bias_1
            )
        )
        # calculate the second residual statistics from the difference between these two biomass estimates (individual for reach triplet and average for all triplets)
        res_2 = forward([agbest_tab_eqbar]) * np.ones((N_lay, 1, 1)) - forward(agbest_tab_eq)
        agbstd_2 = np.nanstd(res_2, axis=(0, 2))
        bias_2 = np.nanmean(res_2, axis=(0, 2))
        # re-calculate biomass estimates averaged over triplets for each individual test using the equation and biases obtained earlier
        agbest_tab_eqbar = inverse(
            (
                np.nansum(curr_output_table[:, :, 4, :], axis=0) / np.nansum(curr_output_table[:, :, 1, :] ** 2, axis=0)
                - bias_1
                - np.array([bias_2]).transpose() * np.ones((1, N_tests))
            )
        )
        # calculate the final biomass estimates by averaging over multiple tests
        agbest = inverse(
            (
                np.nanmean(
                    np.nansum(curr_output_table[:, :, 4, :], axis=0)
                    / np.nansum(curr_output_table[:, :, 1, :] ** 2, axis=0),
                    axis=1,
                )
                - bias_1
                - bias_2
            )
        )
        # calculate the third residual statistics from the difference between the individual test estimates and the final biomass values
        res_3 = forward([agbest]).transpose() * np.ones((1, N_tests)) - forward(agbest_tab_eqbar)
        agbstd_3 = np.nanstd(res_3, axis=(1))
        bias_3 = np.nanmean(res_3, axis=(1))
        # re-calculate the final biomass estimates by averaging over multiple tests and correcting for biases

        agbest = inverse(
            (
                np.nanmean(
                    np.nansum(curr_output_table[:, :, 4, :], axis=0)
                    / np.nansum(curr_output_table[:, :, 1, :] ** 2, axis=0),
                    axis=1,
                )
                - bias_1
                - bias_2
                - bias_3
            )
        )

        # consolidate the final error values

        ## ESTIMATING FINAL PARAMETERS
        # agb estimate table matching the shape of sigma table
        agbest_tab = np.array([agbest]) * np.ones((N_laypol, 1))
        # mask for valid data
        curr_mask = (acqid_tab != -1) & ~np.isnan(agbest_tab)
        # extract observables
        s = forward(sigma_tab[curr_mask])
        c = forward(np.cos(theta_tab[curr_mask]))
        w = forward(np.maximum(1, agbest_tab))[curr_mask]
        # create parameter index vectors
        i_l, i_a, i_n = [par_tab[curr_mask] + par_offset for par_tab, par_offset in zip(par_tabs, par_offsets)]
        # number of parameters to be estimated
        N_par = np.max(i_n) + 1
        # number of observables
        N = len(s)
        # cost function
        def F(x):
            # calculate parameters from x
            l = l0 + dl * np.sin(x[i_l]) ** 2
            a = a0 + da * np.sin(x[i_a]) ** 2
            n = n0 + dn * np.sin(x[i_n]) ** 2
            # calculate cost function
            A = (l + a * w + n * c - s) ** 2
            return 1 / N * np.nansum(A)

        # pre-calculate selection matrices
        X = [
            np.column_stack(([np.float32(i_l == ii) for ii in np.arange(N_par)])),
            np.column_stack(([np.float32(i_a == ii) for ii in np.arange(N_par)])),
            np.column_stack(([np.float32(i_n == ii) for ii in np.arange(N_par)])),
        ]
        # gradient
        def G(x):
            # calculate parameters from x
            l = l0 + dl * np.sin(x[i_l]) ** 2
            a = a0 + da * np.sin(x[i_a]) ** 2
            n = n0 + dn * np.sin(x[i_n]) ** 2
            # calculate parameter derivatives wrt x
            D_l = dl * np.sin(2 * x[i_l])
            D_a = da * np.sin(2 * x[i_a])
            D_n = dn * np.sin(2 * x[i_n])
            # calculate cost function gradient
            A = l + a * w + n * c - s
            return 2 * (
                1
                / N
                * np.nansum(
                    np.array([A * D_l]).transpose() * X[0]
                    + np.array([A * D_a * w]).transpose() * X[1]
                    + np.array([A * D_n * c]).transpose() * X[2],
                    axis=0,
                )
            )

        # initial x-parameter vector (random)
        x0 = np.pi / 2 * np.random.rand(N_par)
        # perform inversion
        mod = sp.optimize.minimize(F, x0, jac=G)
        # show message
        print('\tFinal parameter estimation finished (%s)' % (mod.message))
        # extract parameters and remove unchanged ones
        x = mod.x
        x[np.abs(x - x0) < 1e-4] = np.nan
        # get final parameter values
        parest = [
            l0 + dl * np.sin(x[: par_offsets[1]]) ** 2,
            a0 + da * np.sin(x[par_offsets[1] : par_offsets[2]]) ** 2,
            n0 + dn * np.sin(x[par_offsets[2] : par_offsets[3]]) ** 2,
        ]

        # create parameter look-up table
        parameter_lut = np.column_stack(
            (
                np.array(
                    [
                        par_block_index_curr,
                        EE_par_coordinates[par_block_index_curr],
                        NN_par_coordinates[par_block_index_curr],
                        EE_par_coordinates[par_block_index_curr] + wE_par,
                        NN_par_coordinates[par_block_index_curr] + wN_par,
                    ]
                )
                * np.ones((np.sum(valid_acqs), 1)),
                acq_info_tab[valid_acqs, 0:2],
                np.column_stack([parest[ipar][par_id_luts[ipar]] for ipar in np.arange(3)]),
            )
        )
        # save parameter

        with open(par_path, 'at') as file_par:
            for curr_row in parameter_lut:
                file_par.write(
                    '\n%d\t%.2f\t%.2f\t%.2f\t%.2f\t%d\t%d\t%.2f\t%.2f\t%.2f'
                    % tuple([curr_val for curr_val in curr_row])
                )
            file_par.write('\n')
            file_par.close()

        ## SAVING LOW RESOLUTION DATA: new cals
        # roi resolution values (estimate and stds)
        # for the reshape see also how the rois are meaned in mean_on_rois ()
        roi_vals = [
            agbest.reshape(len(NN_roi_axis), len(EE_roi_axis)),
            (agbest * 0 + agbstd_1).reshape(len(NN_roi_axis), len(EE_roi_axis)),
            (agbest * 0 + agbstd_2).reshape(len(NN_roi_axis), len(EE_roi_axis)),
            (agbest * 0 + agbstd_3).reshape(len(NN_roi_axis), len(EE_roi_axis)),
        ]

        upper_left_easting_coord = EE_roi_axis[0]  # i.e. horizontal
        sampling_step_east_west = EE_roi_axis[1] - EE_roi_axis[0]
        upper_left_northing_coord = NN_roi_axis[0]  # i.e. vertical
        sampling_step_north_south = NN_roi_axis[1] - NN_roi_axis[0]
        geotransform = [
            upper_left_easting_coord,
            sampling_step_east_west,
            0,
            upper_left_northing_coord,
            0,
            sampling_step_north_south,
        ]

        lut_cal_curr_row = [EE_roi_axis[0], EE_roi_axis[-1], NN_roi_axis[-1], NN_roi_axis[0], 1]
        lut_cal = np.vstack([lut_cal, lut_cal_curr_row])

        cal_out_intermediate = os.path.join(
            temp_output_folder_agb_est, 'cal_estimated_from_block_{}'.format(par_block_index_curr)
        )

        curr_cal_path = tiff_formatter(
            roi_vals, cal_out_intermediate, geotransform, multi_layers_tiff=True, gdal_data_format=gdal.GDT_Float32
        )

        lut_cal_paths.append(curr_cal_path)

        sigma_tab = np.load(
            os.path.join(sigma_tab_pixel_folder, 'sigma_tab_pixels_block_{}.npy'.format(par_block_index_curr))
        )
        theta_tab = np.load(
            os.path.join(theta_tab_pixel_folder, 'theta_tab_pixels_block_{}.npy'.format(par_block_index_curr))
        )
        acqid_tab = np.load(
            os.path.join(acq_tab_pixel_folder, 'acq_tab_pixels_block_{}.npy'.format(par_block_index_curr))
        )
        polid_tab = np.load(
            os.path.join(polid_tab_pixel_folder, 'polid_tab_pixels_block_{}.npy'.format(par_block_index_curr))
        )

        par_tabs = [
            tableLookupFloat(
                parameter_lut[:, 5:7],
                parameter_lut[:, 7 + ipar],
                np.column_stack((polid_tab.flatten(), acqid_tab.flatten())),
            ).reshape(polid_tab.shape)
            for ipar in np.arange(3)
        ]
        agbest_tab = inverse(
            (
                np.nansum(
                    (forward(sigma_tab) - par_tabs[0] - par_tabs[2] * forward(np.cos(theta_tab))) * par_tabs[1], axis=0
                )
                / np.nansum(par_tabs[1] ** 2, axis=0)
            )
        )

        # full resolution values (estimate and stds)
        agbstd_2 = agbstd_2.reshape(len(NN_roi_axis), len(EE_roi_axis))
        agbstd_3 = agbstd_3.reshape(len(NN_roi_axis), len(EE_roi_axis))

        pix_vals = [
            agbest_tab.reshape(len(NN_pixel_axis), len(EE_pixel_axis)),
            (agbest_tab * 0 + agbstd_1).reshape(len(NN_pixel_axis), len(EE_pixel_axis)),
            sp.interpolate.griddata(
                (EE_roi_mesh[~np.isnan(agbstd_2)], NN_roi_mesh[~np.isnan(agbstd_2)]),
                agbstd_2[~np.isnan(agbstd_2)],
                (EE_pixel_mesh, NN_pixel_mesh),
                method='nearest',
            )
            * (~np.isnan(agbest_tab.reshape(len(NN_pixel_axis), len(EE_pixel_axis)))),
            sp.interpolate.griddata(
                (EE_roi_mesh[~np.isnan(agbstd_3)], NN_roi_mesh[~np.isnan(agbstd_3)]),
                agbstd_3[~np.isnan(agbstd_3)],
                (EE_pixel_mesh, NN_pixel_mesh),
                method='nearest',
            )
            * (~np.isnan(agbest_tab.reshape(len(NN_pixel_axis), len(EE_pixel_axis)))),
        ]

        quality_layer = np.sqrt(pix_vals[1] ** 2 + pix_vals[2] ** 2 + pix_vals[3] ** 2)

        agb_est_ground_out_file_curr = os.path.join(
            temp_output_folder_agb_est, 'agb_est_block_{}.tif'.format(par_block_index_curr)
        )
        geotransform_curr = [EE_pixel_axis[0], dE, 0, NN_pixel_axis[0], 0, dN]

        agb_est_ground_out_file_curr = tiff_formatter(
            [pix_vals[0], quality_layer],
            agb_est_ground_out_file_curr,
            geotransform_curr,
            gdal_data_format=gdal.GDT_Float32,
            projection=projection_curr,
            multi_layers_tiff=True,
        )

        # [(left, lower), (right, upper)]
        try:
            lon_min, lat_min = getattr(e7g_intermediate, subgrid_code).xy2lonlat(min(EE_pixel_axis), min(NN_pixel_axis))
            lon_max, lat_max = getattr(e7g_intermediate, subgrid_code).xy2lonlat(max(EE_pixel_axis), max(NN_pixel_axis))
        except Exception as e:
            logging.error(
                'Cannot recognize input FNF Equi7 mask "{}" sub-grid folder name :'.format(subgrid_code) + str(e),
                exc_info=True,
            )
            raise

        bbox = [(lon_min, lat_min), (lon_max, lat_max)]
        ftiles = e7g_product.search_tiles_in_roi(bbox=bbox)

        equi7_agb_est_outdir_curr = os.path.join(
            temp_output_folder_agb_est, 'eq7_agb_est_block_{}'.format(par_block_index_curr)
        )
        os.makedirs(equi7_agb_est_outdir_curr)
        equi7_agb_est_outdir_temp = image2equi7grid(
            e7g_product,
            agb_est_ground_out_file_curr,
            equi7_agb_est_outdir_curr,
            gdal_path=gdal_path,
            ftiles=ftiles,
            accurate_boundary=False,
            tile_nodata=np.nan,
        )

        for idx, equi7_tiff_name in enumerate(equi7_agb_est_outdir_temp):

            driver = gdal.Open(equi7_tiff_name, GA_ReadOnly)
            data = driver.ReadAsArray()
            driver = None

            if np.sum(np.isnan(data)) == data.shape[0] * data.shape[1]:
                shutil.rmtree(os.path.dirname(equi7_tiff_name))
            else:
                equi7_agb_est_out_tif_names_not_merged.append(equi7_tiff_name)

        ## PREPARING NEXT TILE
        # swap flag
        par_finished[par_block_index_curr] = True
        # remove current par block from list
        block_order = block_order[block_order != par_block_index_curr]
        # select next par block
        if len(block_order) > 0:
            # if there is at least one left, just take the next closest to CALdata
            par_block_index_curr = block_order[0]
        else:
            break

    ####################### FINAL  EQUI7 TILES MERGING ########################
    logging.info('AGB: final step, merging equi7 blocks togeter...')
    equi7_agb_est_out_tif_names_per_tile = {}
    for equi7_agb_est_out_tif_name in equi7_agb_est_out_tif_names_not_merged:

        tile_name = os.path.basename(os.path.dirname(equi7_agb_est_out_tif_name))

        if tile_name in equi7_agb_est_out_tif_names_per_tile.keys():
            equi7_agb_est_out_tif_names_per_tile[tile_name].append(equi7_agb_est_out_tif_name)
        else:
            equi7_agb_est_out_tif_names_per_tile[tile_name] = [equi7_agb_est_out_tif_name]

    for (tile_name, equi7_agb_est_out_tif_names_curr_tile) in equi7_agb_est_out_tif_names_per_tile.items():

        driver = gdal.Open(equi7_agb_est_out_tif_names_curr_tile[0], GA_ReadOnly)
        data_merged = np.nan * np.zeros((driver.RasterYSize, driver.RasterXSize))
        quality_merged = np.nan * np.zeros((driver.RasterYSize, driver.RasterXSize))
        driver = None

        for idx_tile, equi7_agb_est_out_tif_name_curr_tile in enumerate(equi7_agb_est_out_tif_names_curr_tile):

            # merging current tile blocks togeter:
            driver = gdal.Open(equi7_agb_est_out_tif_name_curr_tile, GA_ReadOnly)

            if idx_tile == 0:
                geotransform_out = driver.GetGeoTransform()
            elif geotransform_out != driver.GetGeoTransform():
                err_str = 'Same equi7 tiles cannot have different geotrasform'
                logging.error(err_str)
                raise ValueError(err_str)

            data_merged = np.nanmean(np.dstack((data_merged, driver.GetRasterBand(1).ReadAsArray())), axis=2)
            quality_merged = np.sqrt(
                np.nanmean(np.dstack((quality_merged ** 2, driver.GetRasterBand(2).ReadAsArray() ** 2)), axis=2)
            )
            driver = None

        invalid_values_mask = np.logical_or(
            data_merged < proc_conf.AGB.estimation_valid_values_limits[0],
            data_merged > proc_conf.AGB.estimation_valid_values_limits[-1],
        )
        data_merged[invalid_values_mask] = np.nan
        quality_merged[invalid_values_mask] = np.nan

        sub_grid_folder_name = os.path.basename(
            os.path.dirname(os.path.dirname(equi7_agb_est_out_tif_names_curr_tile[0]))
        )
        output_folder_curr_tile = os.path.join(
            output_folder, sub_grid_folder_name, tile_name, 'agb_est_' + tile_name + sub_grid_folder_name + '.tif'
        )
        if not os.path.exists(os.path.dirname(output_folder_curr_tile)):
            os.makedirs(os.path.dirname(output_folder_curr_tile))
        output_folder_curr_tile = tiff_formatter(
            [data_merged, quality_merged],
            output_folder_curr_tile,
            geotransform_out,
            gdal_data_format=gdal.GDT_Float32,
            projection=projection_curr,
            multi_layers_tiff=True,
        )
        logging.info('    equi7 tile blocks merged and savad to file{}'.format(output_folder_curr_tile))

    if proc_conf.delete_temporary_files:
        try:
            shutil.rmtree(temp_output_folder)
        except:
            pass
    logging.info('AGB: estimation ended correctly.\n')
