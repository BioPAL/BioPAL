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
from biopal.agb.processing_AGB import (
    LookupTableAGB,
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
from biopal.data_operations.data_operations import (
    read_and_oversample_data,
    read_and_oversample_aux_data,
    fnf_equi7_load_filter_equi7format,
    fnf_tandemx_load_filter_equi7format,
    apply_dem_flattening,
)
from biopal.utility.utility_functions import (
    Task,
    choose_equi7_sampling,
    check_if_path_exists,
    check_if_geometry_auxiliaries_are_present,
    check_fnf_folder_format,
    check_cal_format,
    get_min_time_stamp_repository,
    get_raster_cal_names,
    get_foss_cal_names,
    resolution_heading_correction,
    decode_unique_acquisition_id_string,
    save_breakpoints,
)
from biopal.geocoding.geocoding import (
    geocoding,
    geocoding_init,
)
from biopal.io.xml_io import (
    parse_chains_input_file,
    parse_chains_configuration_file,
)
from biopal.io.data_io import (
    tiff_formatter,
)
from biopal.screen_calibration.screen_calibration import apply_calibration_screens
from biopal.geometry.utility_geometry import compute_and_oversample_geometry_auxiliaries
from biopal.ground_cancellation.ground_cancellation import ground_cancellation


class AboveGroundBiomass(Task):
    def __init__(
        self,
        configuration_file_xml,
        geographic_boundaries,
        geographic_boundaries_per_stack,
        gdal_path,
    ):
        super().__init__(configuration_file_xml)
        self.geographic_boundaries = geographic_boundaries
        self.geographic_boundaries_per_stack = geographic_boundaries_per_stack
        self.gdal_path = gdal_path

    def _run(self, input_file_xml):

        # Main APP #1: Stack Based Processing
        stack_based_processing_obj = StackBasedProcessing(
            self.configuration_file_xml,
            self.geographic_boundaries,
            self.geographic_boundaries_per_stack,
            self.gdal_path,
        )

        # Run Main APP #1: Stack Based Processing
        (lut_cal, lut_fnf, lut_stacks, equi7_initialization) = stack_based_processing_obj.run(
            input_file_xml
        )

        # Main APP #2: AGB Core Processing
        agb_processing_obj = AGBCoreProcessing(
            self.configuration_file_xml,
            self.geographic_boundaries,
            self.gdal_path,
            lut_cal,
            lut_fnf,
            lut_stacks,
            equi7_initialization,
        )

        # Run Main APP #2: AGB Core Processing
        agb_processing_obj.run(input_file_xml)


class StackBasedProcessing(Task):
    def __init__(
        self,
        configuration_file_xml,
        geographic_boundaries,
        geographic_boundaries_per_stack,
        gdal_path,
    ):
        super().__init__(configuration_file_xml)
        self.geographic_boundaries = geographic_boundaries
        self.geographic_boundaries_per_stack = geographic_boundaries_per_stack
        self.gdal_path = gdal_path

    def _run(self, input_file_xml):

        ########################## INITIAL STEPS ##############################

        logging.info('AGB: Reading chains configuration files')
        check_if_path_exists(self.configuration_file_xml, 'FILE')
        proc_conf = parse_chains_configuration_file(self.configuration_file_xml)
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
        time_tag_mjd_initial = get_min_time_stamp_repository(
            proc_inputs.L1c_repository, proc_inputs.stack_composition
        )

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

                lut_cal[idx_cal, :] = [
                    east_min,
                    east_max,
                    min(north_in, north_out),
                    max(north_in, north_out),
                    flag_cal,
                ]
                lut_cal_paths.append(equi7_cal_fname)

                logging.info('using calibration file: ' + equi7_cal_fname)

        ### initialize the equi7 sampling grid
        equi7_sampling_intermediate = choose_equi7_sampling(
            proc_conf.AGB.intermediate_ground_averaging,
            proc_conf.AGB.intermediate_ground_averaging / 2,
        )
        e7g_intermediate = Equi7Grid(equi7_sampling_intermediate)
        logging.info(
            'EQUI7 Grid sampling used for intermediate products: {}'.format(
                equi7_sampling_intermediate
            )
        )

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
                self.gdal_path,
                self.geographic_boundaries,
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
                self.gdal_path,
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

            lut_fnf[fnf_idx, :] = [
                east_min,
                east_max,
                min(north_in, north_out),
                max(north_in, north_out),
            ]
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
            temp_output_folder_gr = os.path.join(
                temp_output_folder, 'ground_range_geometry', unique_stack_id
            )
            temp_output_folder_e7 = os.path.join(
                temp_output_folder, 'ground_equi7_geometry', unique_stack_id
            )
            os.makedirs(temp_output_folder_gr)
            os.makedirs(temp_output_folder_e7)

            ### load data ( and oversample if requested and if needed )
            try:
                logging.info(
                    'AGB: Data loading for stack ' + unique_stack_id + '; this may take a while:'
                )

                (
                    beta0_calibrated,
                    master_id,
                    raster_info,
                    raster_info_orig,
                ) = read_and_oversample_data(
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
                        logging.info(
                            'AGB: calling geometry library for stack ' + unique_stack_id + '\n'
                        )
                        if geometry_aux_are_present:
                            logging.warning(
                                '    geometry auxiliaries will be overwritten for stack '
                                + unique_stack_id
                            )
                        else:
                            logging.info('\n')
                    else:
                        logging.warning(
                            'AGB: calling geometry library since AuxiliaryProductsFolder "Geometry" is empty or not complete \n'
                        )

                    (
                        _,
                        _,
                        ellipsoid_slope,
                        _,
                        _,
                        _,
                        sar_geometry_master,
                    ) = compute_and_oversample_geometry_auxiliaries(
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
                    off_nadir_angle_rad[master_id] = (
                        off_nadir_angle_rad[master_id] - ellipsoid_slope
                    )
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

                    (
                        _,
                        _,
                        _,
                        _,
                        _,
                        _,
                        _,
                        cal_screens,
                        cal_screens_raster_info,
                        _,
                        _,
                        _,
                    ) = read_and_oversample_aux_data(
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
                logging.error(
                    'AGB: error during auxiliary data computation and/or loading: ' + str(e),
                    exc_info=True,
                )
                raise

            ### Screen calibration (ground steering)
            try:
                if proc_conf.apply_calibration_screen:
                    logging.info('AGB: applying calibration screen...')
                    beta0_calibrated = apply_calibration_screens(
                        beta0_calibrated,
                        raster_info,
                        cal_screens,
                        cal_screens_raster_info,
                        master_id,
                    )
                    logging.info('...done.\n')

                elif proc_conf.DEM_flattening:
                    logging.info('AGB: DEM flattening... ')
                    beta0_calibrated = apply_dem_flattening(
                        beta0_calibrated, kz, reference_height, master_id, raster_info
                    )
                    logging.info('...done.\n')

            except Exception as e:
                logging.error(
                    'AGB: error during screen calibration or DEM flattening.' + str(e),
                    exc_info=True,
                )
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

            windtm_x = np.int(
                np.round(sigma_ground_res_m / raster_info.pixel_spacing_az / 2) * 2 + 1
            )
            windtm_y = np.int(
                np.round(
                    sigma_ground_res_m
                    / (raster_info.pixel_spacing_slant_rg / np.sin(look_angle_rad))
                    / 2
                )
                * 2
                + 1
            )

            sub_factor_x = np.int((windtm_x - 1) / 2)
            sub_factor_y = np.int((windtm_y - 1) / 2)

            logging.info('AGB: multilooking of incidence angle...')
            theta_multi_looked_sr = convolve2d(
                off_nadir_angle_rad[master_id] - slope,
                np.ones((windtm_y, windtm_x)) / windtm_y / windtm_x,
                mode='same',
            )
            theta_multi_looked_sr = theta_multi_looked_sr[::sub_factor_y, ::sub_factor_x]
            logging.info('...done.\n')

            sigma0_sr = {}
            for pol_name in DN_beta0_notched.keys():

                logging.info(
                    'AGB: multilooking of ground notched for polarization {}...'.format(pol_name)
                )
                beta0_notched_multi_looked = convolve2d(
                    np.absolute(DN_beta0_notched[pol_name]) ** 2,
                    np.ones((windtm_y, windtm_x)) / windtm_y / windtm_x,
                    mode='same',
                )

                logging.info('AGB: sigma0 computation for polarization {}...'.format(pol_name))
                sigma0_sr[pol_name] = beta0_notched_multi_looked[
                    ::sub_factor_y, ::sub_factor_x
                ] * np.sin(theta_multi_looked_sr)

                sigma0_sr[pol_name][sigma0_sr[pol_name] < 0] = np.NaN

                logging.info('...done.\n')

            del beta0_notched_multi_looked

            ### saving breakpoints
            if proc_conf.save_breakpoints:
                logging.info(
                    'AGB: saving breakpoints (in slant range geometry) on '
                    + breakpoints_output_folder
                )
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
                    sub_factor_x * raster_info.pixel_spacing_az,
                    sub_factor_y * raster_info.pixel_spacing_slant_rg,
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

                    logging.info(
                        'AGB: geocoding the sigma0 for polarization {}...'.format(pol_name)
                    )
                    sigma0_gr[pol_name] = geocoding(
                        sigma0_sr[pol_name],
                        lon_in,
                        lat_in,
                        lonMeshed_out,
                        latMeshed_out,
                        valid_values_mask,
                    )

                    logging.info('...done.\n')

                del sigma0_sr

                # geocode the theta incidence angle
                logging.info('AGB: Geocoding of incidence angle...')
                theta_multi_looked_gr = geocoding(
                    theta_multi_looked_sr,
                    lon_in,
                    lat_in,
                    lonMeshed_out,
                    latMeshed_out,
                    valid_values_mask,
                )
                logging.info('...done.\n')
                del theta_multi_looked_sr

                logging.info('...done.\n')

            except Exception as e:
                logging.error('AGB: error during geocoding: ' + str(e), exc_info=True)
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
                        sigma0_gr[pol_name],
                        sigma0_ground_fnames[pol_name],
                        geotransform,
                        gdal_data_format=gdal.GDT_Float32,
                    )

                del sigma0_gr

                # geotiff of the theta
                theta_ground_fname = os.path.join(
                    temp_output_folder_gr, 'theta' + '_' + unique_stack_id + '.tif'
                )

                tiff_formatter(
                    theta_multi_looked_gr,
                    theta_ground_fname,
                    geotransform,
                    gdal_data_format=gdal.GDT_Float32,
                )

                del theta_multi_looked_gr
                logging.info('...done.\n')

            except Exception as e:
                logging.error('AGB: error during GEOTIFF formatting: ' + str(e), exc_info=True)
                raise

            ### formatting data to EQUI7
            logging.info(unique_stack_id + ': formatting into EQUI7 grid...')
            try:

                sigma0_equi7_fnames[unique_stack_id] = {}
                # equi7 of the sigma0 (three polarizations)
                for pol_name in sigma0_ground_fnames.keys():

                    equi7_sigma0_outdir = os.path.join(temp_output_folder_e7, 'sigma0_' + pol_name)

                    logging.info(
                        'image2equi7grid IN: '
                        + sigma0_ground_fnames[pol_name]
                        + ' , OUT:'
                        + equi7_sigma0_outdir
                    )
                    sigma0_equi7_fnames[unique_stack_id][pol_name] = image2equi7grid(
                        e7g_intermediate,
                        sigma0_ground_fnames[pol_name],
                        equi7_sigma0_outdir,
                        gdal_path=self.gdal_path,
                        inband=None,
                        subgrid_ids=None,
                        accurate_boundary=False,
                        resampling_type='bilinear',
                        tile_nodata=np.nan,
                    )

                # equi7 of the theta
                equi7_theta_outdir = os.path.join(temp_output_folder_e7, 'theta')

                logging.info(
                    'image2equi7grid IN: ' + theta_ground_fname + ' , OUT:' + equi7_theta_outdir
                )
                theta_equi7_fnames[unique_stack_id] = image2equi7grid(
                    e7g_intermediate,
                    theta_ground_fname,
                    equi7_theta_outdir,
                    gdal_path=self.gdal_path,
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
                raise RuntimeError(error_message + ' :' + theta_equi7_fnames[unique_stack_id][0])

            # N0: upper left corner for output map (UTM 32S)
            _, x_upper_left, y_upper_left = e7g_intermediate.lonlat2xy(
                self.geographic_boundaries_per_stack[unique_stack_id].lon_min,
                self.geographic_boundaries_per_stack[unique_stack_id].lat_max,
            )
            north_in = y_upper_left
            east_min = x_upper_left

            # lower right corner for output map (UTM 32S)
            _, x_lower_right, y_lower_right = e7g_intermediate.lonlat2xy(
                self.geographic_boundaries_per_stack[unique_stack_id].lon_max,
                self.geographic_boundaries_per_stack[unique_stack_id].lat_min,
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

        equi7_initialization = {
            'equi7_sampling_intermediate': equi7_sampling_intermediate,
            'e7g_intermediate': e7g_intermediate,
        }

        ########################## END OF STACK BASED STEPS ######################
        return (
            LookupTableAGB(paths=lut_cal_paths, boundaries=lut_cal, progressive=None),
            LookupTableAGB(paths=lut_fnf_paths, boundaries=lut_fnf, progressive=None),
            LookupTableAGB(
                paths=lut_stacks_paths,
                boundaries=lut_stacks_boundaries,
                progressive=lut_progressive_stacks,
            ),
            equi7_initialization,
        )


class AGBCoreProcessing(Task):
    def __init__(
        self,
        configuration_file_xml,
        geographic_boundaries,
        gdal_path,
        lut_cal,
        lut_fnf,
        lut_stacks,
        equi7_conf=None,
    ):
        if equi7_conf is None:
            equi7_conf = {}
        super().__init__(configuration_file_xml)
        self.geographic_boundaries = geographic_boundaries
        self.gdal_path = gdal_path
        self.lut_cal_paths = lut_cal.paths
        self.lut_cal = lut_cal.boundaries
        self.lut_fnf_paths = lut_fnf.paths
        self.lut_fnf = lut_fnf.boundaries
        self.lut_stacks_paths = lut_stacks.paths
        self.lut_stacks_boundaries = lut_stacks.boundaries
        self.lut_progressive_stacks = lut_stacks.progressive
        self.equi7_sampling_intermediate = equi7_conf.get('equi7_sampling_intermediate')
        self.e7g_intermediate = equi7_conf.get('e7g_intermediate')

    def check_auxiliaries(self, proc_conf):
        if not self.equi7_sampling_intermediate:
            # initialize the equi7 sampling grid
            self.equi7_sampling_intermediate = choose_equi7_sampling(
                proc_conf.AGB.intermediate_ground_averaging,
                proc_conf.AGB.intermediate_ground_averaging / 2,
            )
            if self.equi7_sampling_intermediate:
                logging.warning(
                    'Redefining EQUI7 Grid sampling used for intermediate products: {}'.format(
                        self.equi7_sampling_intermediate
                    )
                )
            else:
                logging.info(
                    'EQUI7 Grid sampling used for intermediate products: {}'.format(
                        self.equi7_sampling_intermediate
                    )
                )
            self.e7g_intermediate = Equi7Grid(self.equi7_sampling_intermediate)
    # %%
    def _run(self, input_file_xml):

        # status update and reading xml files
        logging.info('AGB: Reading chains configuration files')
        check_if_path_exists(self.configuration_file_xml, 'FILE')
        proc_conf = parse_chains_configuration_file(self.configuration_file_xml)
        proc_inputs = parse_chains_input_file(input_file_xml)

        # setting up directories and making sure that preprocessing has been run
        products_folder = os.path.join(proc_inputs.output_folder, 'Products')
        temp_proc_folder = os.path.join(products_folder, 'temporary_processing')
        if not (os.path.exists(temp_proc_folder)):
            error_message = '"temporary_processing" folder is not present in output: StackBasedProcessing APP should be launched before AGBCoreProcessing '
            logging.error(error_message)
            raise RuntimeError(error_message)
        # check auxiliaries (equi7 initialization) and if not present, compute them
        self.check_auxiliaries(proc_conf)

        # preparing folders for final and temporary output
        global_agb_folder = os.path.join(products_folder, 'global_AGB')
        temp_agb_folder = os.path.join(temp_proc_folder, 'agb_estimation')
        os.makedirs(global_agb_folder)
        os.makedirs(temp_agb_folder)
        
        
        ### don't think we need those - it's not reasonable to average with previous estimates
        # sigma_tab_pixel_folder = os.path.join(temp_agb_folder, 'sigma_tab')
        # os.makedirs(sigma_tab_pixel_folder)
        # theta_tab_pixel_folder = os.path.join(temp_agb_folder, 'theta_tab')
        # os.makedirs(theta_tab_pixel_folder)
        # acq_tab_pixel_folder = os.path.join(temp_agb_folder, 'acq_tab')
        # os.makedirs(acq_tab_pixel_folder)
        # polid_tab_pixel_folder = os.path.join(temp_agb_folder, 'polid_tab')
        # os.makedirs(polid_tab_pixel_folder)

        # # parameter text file path
        # par_path = os.path.join(global_agb_folder, 'par_est.txt')

        # read and initialize all the parameters needed for the inversion
        ### future improvement: update the function below and the xml file to incorporate the new structure
        (
            _,
            first_pixel_north,
            first_pixel_east,
            last_pixel_north,
            last_pixel_east,
            pixel_size_east,
            pixel_size_north,
            sample_spacing_east,
            sample_spacing_north,
            _,
            _,
            block_spacing_east,
            block_spacing_north,
            block_size_east,
            block_size_north,
            parameter_variabilities_orig,
            parameter_limits_l,
            parameter_limits_a,
            parameter_limits_n,
            parameter_limits_w,
            number_of_subsets,
            geographic_grid_sampling,
            sub_grid_string,
        ) = initialize_inversion_parameters(
            self.equi7_sampling_intermediate,
            proc_inputs.geographic_grid_sampling,
            self.geographic_boundaries,
            proc_conf.AGB,
        )
        
        ### new functions, to be moved out of this method

        # this function converts the given formula using the list of observable names to
        # a formula using different columns of an array with name in table name
        # note: if multiple columns have the same name (e.g., in the case of agb, where agb is both an observable and a parameter),
        #    the the column position_if_repeated is taken; this way, this function can be used for both calibration data and estimation data
        def convert_formula(in_string,table_name,column_names,position_if_repeated):
            out_string = '0'
            for independent_terms in in_string.split('&'):
                additive_terms_string = '0'
                for additive_terms in independent_terms.split('+'):
                    multiplicative_terms_string = '1'
                    for multiplicative_term in additive_terms.split('*'):
                        curr_position = np.where(multiplicative_term.strip()==np.array(column_names))[0]
                        if len(curr_position)>1:
                            column_number = curr_position[position_if_repeated]
                        else:
                            column_number = curr_position[0]
                        multiplicative_terms_string = '%s * %s[:,%d]' % (multiplicative_terms_string,table_name,column_number)
                    additive_terms_string = '%s + %s' % (additive_terms_string,multiplicative_terms_string) 
                out_string = '%s + (%s)**2' % (out_string,additive_terms_string)
            return out_string.replace('0 + ','').replace('1 * ','')
        
        # define transform functions
        def transform_function(in_data,interval,kind,do_forward=True):
            out_data = np.copy(in_data)
            out_data[(out_data<interval[0]) | (out_data>interval[1])] = np.nan
            # note: no check of data and kind is done here, 
            # it is assumed that the inputs are correct
            if kind=='db':
                if do_forward:
                    out_data = 10*np.log10(out_data)
                else:
                    out_data = 10**(0.1*in_data)
            elif kind=='-db':
                if do_forward:
                    out_data = -10*np.log10(in_data)
                else:
                    out_data = 10**(-0.1*in_data)
            elif kind=='-2db':
                if do_forward:
                    out_data = -10*np.log10(2*in_data)
                else:
                    out_data = 10**(-0.1*in_data)/2
            elif kind=='cosdb':
                if do_forward:
                    out_data = 10*np.log10(np.cos(in_data))
                else:
                    out_data = np.arccos(10**(0.1*in_data))
            else:
                out_data = np.copy(in_data)
            return out_data
        # define parameter transfer functions
        def parameter_transfer_function(x,pmin,pmax,x_is_p=False):
            # note: no check of x, xmax, xmin and kind is done here, 
            # it is assumed that the inputs are correct
            if not x_is_p:
                return pmin+(pmax-pmin)*np.sin(x)**2
            else:
                return np.arcsin(np.sqrt((x-pmin)/(pmax-pmin)))
    
        

        # function for calculating a given statistic for polygons of arbitrary shape  
        def stats_on_cal_polygons(data,pixel_axis_east,pixel_axis_north,cal_polygons,method):
            # here, create a function that calculates statistic in method on CAL data polygons
            #
            print('not implemented')
            return []
        # function for calculating a given statistic on all samples on a grid and all polygons
        def stats_on_all_samples(data,pixel_axis_east,pixel_axis_north,sampling_axis_east,sample_size_east,sampling_axis_north,sample_size_north,cal_polygons,method):
            stats = mean_on_rois(data,pixel_axis_east,pixel_axis_north,sampling_axis_east,sample_size_east,sampling_axis_north,sample_size_north,method)
            if cal_polygons:
                stats_cal_polygons = stats_on_cal_polygons(data,pixel_axis_east,pixel_axis_north,cal_polygons,method)
                stats = np.concatenate((stats,stats_cal_polygons))
            return stats
        
        
        # function for creating list of format strings and a header for subsequent saving of tables into text files
        def get_fmt_and_header(column_names,all_column_groups,all_data_types,delimiter='\t',precision=2,column_width=10):
            out_format = []
            out_header = []
            for curr_column in column_names:
                for curr_column_group,curr_data_type in zip(all_column_groups,all_data_types):
                    if curr_column in curr_column_group:
                        if curr_data_type == 'd':
                            out_format.append('%s%dd' % (r'%',column_width))
                        elif curr_data_type == 'f':
                            out_format.append('%s%d.%df' % (r'%',column_width,precision))
                        break
                out_header.append('%s%ds' % (r'%',column_width) % curr_column)
            return out_format,delimiter.join(out_header)
        
        ########## PREPARING INPUT
        # this should be in the xml file and come out from the new version of the xml read function
        
        
        ### in this proposed new version, the user defines the residual function formula on the following format:
        # formula = 'logsigmahh ~ constanthh + agb_scalinghh * logagb + angle_scalinghh * logcostheta & 
        #               ... &
        #            logsigmavv ~ constantvv + agb_scalingvv * logagb + angle_scalingvv * logcostheta ']
        # Although this format may appear complex at first, it allows a lot of flexibility in defining (multiple) cost functions
        # With this, different models can be used for different polarisations and we can also change between forward and inverse model fitting
        #
        # Four signs are allowed:
        #   1 - &: separation of independent cost functions
        #   2 - ~: equivalent to subtraction of whatever follows
        #   3 - +: addition
        #   4 - *: multiplication
        #   The formula is split based on the order above and the mathematical operations are done in the opposite order
        #     here, we hardcode this before the xml and xml reading function above are updated
        formula = ' & '.join(['l_hh + a_hh * agb_1_db + n_hh * cos_local_db + neg_sigma0_hh_db',
                                'l_hv + a_hv * agb_1_db + n_hv * cos_local_db + neg_sigma0_hv_db',
                                'l_vv + a_vv * agb_1_db + n_vv * cos_local_db + neg_sigma0_vv_db'])
        
            
        ### then, the user has to define each of the elements of the formula
        # Essentially, the unique elements of formula_split are cycled through one by one and the user has to select whether it is an "observable" (read only data), "parameter" (estimate only data),
        # or both ("read and write", typical for AGB for which some calibration data are required)
        # here, we hardcode this thing before the xml file and "initialize_inversion_parameters" are updated
        # first, the names of the observables in formula are included (plus a few extras, which may be needed later, like other agb-related parametrs in db)
        observable_names = ['neg_sigma0_hh_db','neg_sigma0_hv_db','neg_sigma0_vv_db','cos_local_db','agb_1_db','agb_2_db','agb_3_db','agb_4_db']
        # then, the observable sources are selected (this is intimately bound to the available in data, so it is based on the stack-based processing output and calibration data availability)
        observable_sources = ['sigma0_hh','sigma0_vh','sigma0_vv','theta','agb_1','agb_2','agb_3','agb_4']
        # permissible intervals for sources
        observable_ranges = [[1e-8,20],[1e-8,20],[1e-8,20],[20*np.pi/180,60*np.pi/180],[1,700],[0,10],[0,10],[0,10]]
        # then, the applied transforms are selected (see transform_function() above)
        observable_transforms = ['-db','-2db','-db','cosdb','db','none','none','none']
        # finally, averaging methods are selected (for averaging pixels within samples)
        # although oftentimes 'mean' is used, this may not always be the case 
        # e.g., for slope aspect angle v, arctan(mean(sin(v))/mean(cos(v))) is a better averaging method as it is not as susceptible to 2pi ambiguities
        observable_averaging_methods = ['mean','mean','mean','mean','mean','mean','mean','mean'] 
        # then, the names of parameters are specified (should match those in formula)
        # note that agb is special: it is defined as both "observable" and "parameter", so it will be 
        # both read from external data and estimated 
        parameter_names = ['l_hh','l_hv','l_vv','a_hh','a_hv','a_vv','n_hh','n_hv','n_vv','agb_1_db']
        # then, parameter limits are selected; we here simply reorganise the limits in xml to the selected format
        parameter_limits = [parameter_limits_l,parameter_limits_l,parameter_limits_l,
                            parameter_limits_a,parameter_limits_a,parameter_limits_a,
                            parameter_limits_n,parameter_limits_n,parameter_limits_n,
                            parameter_limits_w]
        # parameter variabilities are defined over eight different dimensions:
        #   spatial - between samples
        #   forest class - between forest classes
        #   stack - between stacks, independent from where they come (equivalent to flagging all remaining flags)
        #   global cycle - between stacks if they come from different global cycle (equivalent to temporal of multiple months)
        #   heading - between stacks if they come from different headings
        #   swath - between stacks if they come from different 150-km swaths (equivalent to temporal of several weeks)
        #   subswath - between stacks if they come from different 50-km subswaths (equivalent to temporal of several days)
        #   azimuth - betwen stacks if they come from different 150-km images overlapping in azimuth (equivalent to temporal of a few seconds)
        # as before, some of the new things will be hardcoded for now, but they should come from the xml
        # compared with before, we need to add two extra columns (spatial and forest class variability) and remove the polarisation column,
        parameter_variabilities_orig = np.column_stack((np.ones((3,2))==0,np.array(parameter_variabilities_orig)[:,1:]))
        # moreover, we need to repeat each parameter three times to represent the three polarisations (which now can have different cost functions)
        parameter_variabilities_orig = parameter_variabilities_orig[np.array([0,0,0,1,1,1,2,2,2]),:]
        # then, we convert from array to list
        parameter_variabilities = []
        for parameter_variability in parameter_variabilities_orig:
            parameter_variabilities.append(parameter_variability)
        # finally, we add variability vector for agb1db so that it matches the order from above
        parameter_variabilities.append(np.array([True,False,False,False,False,False,False,False]))
        # derived quantities
        number_of_observables = len(observable_names)
        number_of_parameters = len(parameter_names)
        ###################### END OF PREPARING INPUT 
        
        
        
        # ### checking the formula and list of observables        
        # # split formula to a nested list of dictionaries where different operations are done at different levels
        # calibration_formula_terms = parse_formula(formula,observable_names+parameter_names,True)
        # estimation_formula_terms = parse_formula(formula,observable_names+parameter_names,False)
        # # to be implemented: check that all terms of the formula are either observables or parameters
       
            
        
        # initialize the equi7 sampling grid or output product
        e7g_product = Equi7Grid(geographic_grid_sampling)
        logging.info(
            'EQUI7 Grid sampling used for final products: {}'.format(geographic_grid_sampling)
        )

        equi7_subgrid_names_list = []
        for lut_stacks_path in self.lut_stacks_paths:
            for equi7_subgrid_name_curr in os.listdir(os.path.join(lut_stacks_path, 'sigma0_hh')):

                equi7_subgrid_names_list.append(equi7_subgrid_name_curr)

        if len(set(equi7_subgrid_names_list)) > 1:
            err_str = 'AGB: biomass L2 prototype does not support multiple equi7 sub grid regions'
            logging.error(err_str)
            raise
        subgrid_code = equi7_subgrid_names_list[0][6:8]
        # projection_prev = ''
        equi7_agb_est_out_tif_names_not_merged = []

       


        # parameter block mesh and flattened coordinate vectors
        block_mesh_east, block_mesh_north = np.meshgrid(np.arange(first_pixel_east, last_pixel_east, block_spacing_east), np.arange(first_pixel_north, last_pixel_north, block_spacing_north))
        block_corner_coordinates_east, block_corner_coordinates_north = [
            np.float64(x.flatten()) for x in [block_mesh_east, block_mesh_north]
        ]
        # count the blocks and set up a "finished" flag vector
        number_of_blocks = len(block_corner_coordinates_east.flatten())
        block_finished_flag = np.zeros((number_of_blocks), dtype='bool')
        
        # read acquisition info table
        stack_info_table = self.lut_progressive_stacks
        stack_info_table_columns = ['stack_id','global_cycle_id','heading_degrees','swath_id','subswath_id','azimuth_id']
        stack_id_vector = stack_info_table[:,0]
        number_of_stacks = len(stack_id_vector)
        
        
        # Compute the starting block and the nblock processing order
        # proposed improvement: use SAR image boundaries in this
        current_block_index, block_order = compute_processing_blocs_order(
            self.lut_cal, block_corner_coordinates_east, block_size_east, block_corner_coordinates_north, block_size_north
        )

        
        # run as long as not all blocks are finished
        counter_blocks_run = 0
        while ~np.all(block_finished_flag):

            # show status
            counter_blocks_run = counter_blocks_run + 1
            logging.info('AGB: Running block {} out of {} (block ID: {})'.format(counter_blocks_run, number_of_blocks, current_block_index))
            
            
            
            
            # extent of current block # [min_east, max_east, min_north, max_north]
            current_block_extents = np.array(
                [
                    block_corner_coordinates_east[current_block_index],
                    block_corner_coordinates_east[current_block_index] + block_size_east,
                    block_corner_coordinates_north[current_block_index],
                    block_corner_coordinates_north[current_block_index] + block_size_north,
                ]
            )  

            # sampling axes and meshes for current block
            sampling_axis_east = np.arange(current_block_extents[0], current_block_extents[1], sample_spacing_east)
            sampling_axis_north = np.arange(current_block_extents[2], current_block_extents[3], sample_spacing_north)
            sampling_mesh_east, sampling_mesh_north = np.meshgrid( sampling_axis_east,
                                                                  sampling_axis_north)
            number_of_samples_on_grid = len(sampling_mesh_east.flatten())
            
            ### to be implemented:
            # define or/and load additional, arbitrarily shaped sampling polygons
            additional_sampling_polygons = []
            number_of_polygon_samples = len(additional_sampling_polygons)
            
            # calculate the total number of samples
            number_of_samples = number_of_samples_on_grid + number_of_polygon_samples
            
            # checking the number of samples
            if number_of_samples_on_grid < proc_conf.AGB.min_number_of_rois:
                logging.info(
                    'skipping block #{} because the number of samples #{} cannot be less than #{}'.format(
                        current_block_index, number_of_samples_on_grid, proc_conf.AGB.min_number_of_rois
                    )
                )
                # swap the flag
                block_finished_flag[current_block_index] = True
                # remove current block from list
                block_order = block_order[block_order != current_block_index]
                # select next block
                if len(block_order) > 0:
                    # if there is at least one left, just take the next closest to CALdata
                    current_block_index = block_order[0]
                    continue
                else:
                    break

            # pixel axes for current block
            pixel_axis_east = np.arange(current_block_extents[0], current_block_extents[1], pixel_size_east)
            pixel_axis_north = np.arange(current_block_extents[2], current_block_extents[3], pixel_size_north)
            pixel_mesh_east, pixel_mesh_north = np.meshgrid( pixel_axis_east, 
                                                            pixel_axis_north)
            number_of_pixels = len(pixel_mesh_east.flatten())

            
            ## THE FOLLOWING PART IS TABULATING
            # could be moved to a separate function
            try:
                
                

                # cycle through stacks and check that there are some data within current block
                block_has_data = np.zeros(number_of_stacks, dtype='bool')
                for stack_idx, stack_path in enumerate(self.lut_stacks_paths):

                    # go ahead only if current parameter block is at least partially contained in the data stack:
                    block_has_data[stack_idx] = check_intersection(
                        self.lut_stacks_boundaries[stack_idx, 0],
                        self.lut_stacks_boundaries[stack_idx, 1],
                        self.lut_stacks_boundaries[stack_idx, 2],
                        self.lut_stacks_boundaries[stack_idx, 3],
                        current_block_extents[0],
                        current_block_extents[1],
                        current_block_extents[3],
                        current_block_extents[2],
                    )
                number_of_stacks_inside = np.sum(block_has_data == True)
                if (number_of_stacks_inside == 0):
                    logging.info(
                        'skipping block #{} due to lack of valid data points'.format(
                            current_block_index
                        )
                    )
                    # swap the flag
                    block_finished_flag[current_block_index] = True
                    # remove current par block from list
                    block_order = block_order[block_order != current_block_index]
                    # select next par block
                    if len(block_order) > 0:
                        # if there is at least one left, just take the next closest to CALdata
                        current_block_index = block_order[0]
                        continue
                    else:
                        break

                # cycle through cals and see that there are some cals within current block
                block_has_cal = np.zeros(len(self.lut_cal_paths), dtype='bool')
                for cal_idx, cal_path in enumerate(self.lut_cal_paths):

                    # go ahead only if current parameter block is at least partially contained in the data stack:
                    block_has_cal[cal_idx] = check_intersection(
                        self.lut_cal[cal_idx, 0],
                        self.lut_cal[cal_idx, 1],
                        self.lut_cal[cal_idx, 2],
                        self.lut_cal[cal_idx, 3],
                        current_block_extents[0],
                        current_block_extents[1],
                        current_block_extents[3],
                        current_block_extents[2],
                    )

                if not np.any(block_has_cal):
                    logging.info(
                        'skipping block #{} due to lack of calibration data'.format(
                            current_block_index
                        )
                    )
                    # swap the flag
                    block_finished_flag[current_block_index] = True
                    # remove current par block from list
                    block_order = block_order[block_order != current_block_index]
                    # select next par block
                    if len(block_order) > 0:
                        # if there is at least one left, just take the next closest to CALdata
                        current_block_index = block_order[0]
                        continue
                    else:
                        break
                # %%    
                    
                    
                # defining names for identifiers (sampleID & forest class ID and then all columns from stack info table)
                identifier_names = ['sample_id','forest_class_id'] + stack_info_table_columns
                number_of_identifiers = len(identifier_names)
                identifier_table = np.nan * np.zeros (( 
                    number_of_samples * number_of_stacks,
                    number_of_identifiers
                    ))
                
                
                    
                # filling out the first column with sample IDs
                identifier_table[:,0] = np.kron(np.arange(number_of_samples),np.ones(number_of_stacks))
                
                # filling out columns 3-8 with stack IDs and corresponding other identifications from stack_info_table
                identifier_table[:,2] = np.kron(np.ones(number_of_samples),stack_id_vector)
                for id_idx in range(5):
                    identifier_table[:,3+id_idx] = sp.interpolate.interp1d(stack_id_vector,stack_info_table[:,1+id_idx],kind='nearest')(identifier_table[:,2])
                
                
                ### READING AND SAMPLING FOREST CLASS DATA
                logging.info('reading & sampling forest class data')
                # get equi7 information
                equi7_grid_name = os.listdir(os.path.join(self.lut_stacks_paths[0], 'sigma0_hh'))[0]
                # loading all forest class maps that fall inside a block: note, more than one equi7 tile can fall inside
                # loaded masks needs to be re-interpolated over the pixels grid (pixel_axis_east, pixel_axis_north)
                forest_class_map_interp = np.zeros(
                    (len(pixel_axis_north), len(pixel_axis_east)), dtype='float'
                )

                counter_forest_class_maps = 0
                for forest_class_map_idx, forest_class_map_path in enumerate(self.lut_fnf_paths):

                    forest_class_map_boundaries = self.lut_fnf[forest_class_map_idx]

                    forest_class_map_is_inside = check_intersection(
                        forest_class_map_boundaries[0],
                        forest_class_map_boundaries[1],
                        forest_class_map_boundaries[2],
                        forest_class_map_boundaries[3],
                        current_block_extents[0],
                        current_block_extents[1],
                        current_block_extents[3],
                        current_block_extents[2],
                    )
                    if forest_class_map_is_inside:

                        forest_class_map_interp_curr = np.round(
                            interp2d_wrapper(
                                forest_class_map_path, 1, pixel_axis_east, pixel_axis_north, fill_value=float(0)
                            )
                        )

                        # mean all the fnf tiles
                        forest_class_map_interp = np.ceil(
                            merge_agb_intermediate(
                                forest_class_map_interp, forest_class_map_interp_curr, method='nan_mean'
                            )
                        )
                        
                        # set all unrealistic values to 0 = non-forest
                        forest_class_map_interp[(forest_class_map_interp <= 0) | np.isnan(forest_class_map_interp)] = 0
                        
                        counter_forest_class_maps = counter_forest_class_maps + 1

                if counter_forest_class_maps == 0:

                    err_str = 'Cannot find any forest class map falling in current block coordinates.'
                    logging.error(err_str)
                    raise ImportError(err_str)
                
                # sampling forest class map
                temp_forest_class_vector = np.int32(np.round(stats_on_all_samples(
                    forest_class_map_interp,
                    pixel_mesh_east,
                    pixel_mesh_north,
                    sampling_axis_east,
                    sample_spacing_east,
                    sampling_axis_north,
                    sample_spacing_north,
                    additional_sampling_polygons,
                    'mode',
                )))
                # repeating it across stacks and inserting into observable data table (note: forest class assumed constant across stacks)
                identifier_table[:,1] = np.kron(temp_forest_class_vector,np.ones(number_of_stacks))
                
                
                
                
                
                
                ### READING OBSERVABLE DATA
                # allocate observable table
                observable_table = np.nan * np.zeros((
                    number_of_samples * number_of_stacks,
                    number_of_observables))
                # the following assumes a certain format on the radar observables (which consist of multiple stacks)
                # and on agb (which is assumed to be the same for all stacks)
                # in case we want to work with multiple agb measurements over time (e.g., agb change), then the following needs to be modified
                observable_types = np.nan * np.zeros(number_of_observables)
                # cycle over stacks
                for stack_idx, stack_path in enumerate(self.lut_stacks_paths):

                    # go ahead only if current stack is (at least partially) contained in the current parameter block:
                    if block_has_data[stack_idx]:

                        # cycle through observables sources in the stack
                        for source_idx,source_name in enumerate(observable_sources):
                            
                            # only do the following for the radar observables (i.e. not an agb observable)
                            if not (source_name.split('_')[0]=='agb'):
                                
                                observable_types[source_idx] = 0
                                
                                source_eq7_parent_dir_name = os.path.join(
                                    stack_path, source_name, equi7_grid_name
                                )
                                equi7_tiles_names = os.listdir(source_eq7_parent_dir_name)
    
                                # cycle over all the equi7 tiles, interpolate over pixel grid and average
                                source_data_interp = np.NaN * np.zeros(
                                    (len(pixel_axis_north), len(pixel_axis_east)), dtype='float'
                                )
                                for equi7_tile_name in equi7_tiles_names:
    
                                    source_tiff_name = os.path.join(
                                        source_eq7_parent_dir_name,
                                        equi7_tile_name,
                                        source_name 
                                        + '_'
                                        + os.path.basename(stack_path)
                                        + '_'
                                        + equi7_grid_name[6:]
                                        + '_'
                                        + equi7_tile_name
                                        + '.tif',
                                    )
    
                                    source_data_interp_curr = interp2d_wrapper(
                                        source_tiff_name,
                                        1,
                                        pixel_axis_east,
                                        pixel_axis_north,
                                        fill_value=np.NaN,
                                    )
    
                                    source_data_interp = merge_agb_intermediate(
                                        source_data_interp, source_data_interp_curr, method='nan_mean'
                                    )
    
                                # masking the stack:
                                source_data_interp[forest_class_map_interp == 0] = np.NaN
        
                                logging.info("sampling and transforming source '{}' to observable '{}' using transform '{}', stack {}, and averaging method '{}' (mean value: {})".format(
                                    source_name,
                                    observable_names[source_idx],
                                    observable_transforms[source_idx],
                                    stack_id_vector[stack_idx],
                                    observable_averaging_methods[source_idx],
                                    np.nanmean(source_data_interp)))
                                
                                # calculate sample statistics
                                temp_transformed_sampled_data = transform_function(
                                    stats_on_all_samples(
                                        source_data_interp,
                                        pixel_mesh_east,
                                        pixel_mesh_north,
                                        sampling_axis_east,
                                        sample_spacing_east,
                                        sampling_axis_north,
                                        sample_spacing_north,
                                        additional_sampling_polygons, 
                                        observable_averaging_methods[source_idx],
                                        ),
                                    observable_ranges[source_idx],
                                    observable_transforms[source_idx])
                                # find rows where the stack id matches the current stack id
                                current_rows = (identifier_table[:,2]==stack_id_vector[stack_idx])
                                # fill out the table
                                observable_table[current_rows,
                                                      source_idx] = temp_transformed_sampled_data
                # go through observable sources again
                for source_idx,source_name in enumerate(observable_sources):
                    
                    # do the following only for agb-type observables
                    if source_name.split('_')[0]=='agb':
                        
                        observable_types[source_idx] = 1
                        
                        agb_data_interp = np.nan * np.zeros(
                            (len(pixel_axis_north), len(pixel_axis_east)), dtype='float'
                        )
                        
                        agb_source_band = int(source_name.split('_')[1])
                    
                        for cal_idx, cal_path in enumerate(self.lut_cal_paths):
                            if block_has_cal[cal_idx]:
                                # counter_cal = counter_cal + 1
                                # logging.info('reading cal raster from path {}'.format(cal_path))
                            
                                agb_data_interp_curr = interp2d_wrapper(
                                    cal_path, agb_source_band, pixel_axis_east, pixel_axis_north, fill_value=np.NaN
                                )
                                # mean all the equi7 tiles
                                agb_data_interp = merge_agb_intermediate(
                                    agb_data_interp, agb_data_interp_curr, method='nan_mean'
                                )
                        logging.info("sampling and transforming source '{}' to observable '{}' using transform '{}', and averaging method '{}' (mean value: {})".format(
                                    source_name,
                                    observable_names[source_idx],
                                    observable_transforms[source_idx],
                                    observable_averaging_methods[source_idx],
                                    np.nanmean(agb_data_interp)))
                        # statistics over samples
                        temp_transformed_sampled_data = transform_function(
                            stats_on_all_samples(
                                agb_data_interp,
                                pixel_mesh_east,
                                pixel_mesh_north,
                                sampling_axis_east,
                                sample_spacing_east,
                                sampling_axis_north,
                                sample_spacing_north,
                                additional_sampling_polygons,
                                method='nan_mean',
                            ),
                            observable_ranges[source_idx],
                            observable_transforms[source_idx])
                        # fill out the table
                        observable_table[:,source_idx] = np.kron(temp_transformed_sampled_data,np.ones(number_of_stacks))
                        
                    
                # for simpler management, create vectors indicating columns with identifiers, observables, and parameter positions, etc
                columns_with_sar_observables = np.where(observable_types==0)[0]
                columns_with_agb_observables = np.where(observable_types==1)[0]
                
                # mark rows in observable data table that have negative identifiers, nan-valued sar observables, infinite sar observables, or negative agb values
                invalid_rows = np.any(identifier_table<0,axis=1) | \
                    np.any(np.isnan(observable_table[:,columns_with_sar_observables]),axis=1) | \
                    np.any(~np.isfinite(observable_table[:,columns_with_sar_observables]),axis=1) | \
                        np.any(observable_table[:,columns_with_agb_observables]<0,axis=1)
                # exclude invalid rows
                observable_table = observable_table[~invalid_rows,:]
                identifier_table = identifier_table[~invalid_rows,:]
                # number of rows in data table
                number_of_rows_in_observable_table = observable_table.shape[0]
                        
                        
                ### PREPARING PARAMETER TABLES
                parameter_property_names = ['lower_limit','upper_limit','initial_value']+['estimate_%d' % (ii) for ii in np.arange(number_of_subsets)+1]
                parameter_position_names = ['row_'+parameter_name for parameter_name in parameter_names]
                parameter_tables = []
                parameter_table_columns = []
                parameter_position_table = np.nan * np.zeros((
                    number_of_rows_in_observable_table,number_of_parameters))
                # creating parameter matrices
                for parameter_idx,parameter_variability in enumerate(parameter_variabilities):
                    # take out only the relevant identifiers (the ones that change as per parameter variability)
                    # and create column names by adding four additional columns: min, max and initial value, and estimated value (later set to NaN)
                    parameter_table_columns.append(np.concatenate((np.array(identifier_names)[parameter_variability],
                                                                   np.array(parameter_property_names))))
                    # create the minimal ID table (without unnecessary columns for those dimension across which the parameter doesn't change)
                    temp_ids_table = identifier_table[:,np.where(parameter_variability)[0]]
                    # create the last four columns
                    temp_minmax_table = np.array([parameter_limits[parameter_idx]]) * np.ones((number_of_rows_in_observable_table,1))
                    temp_initial_table = np.mean(temp_minmax_table,axis=1)
                    temp_estimated_table = np.kron(np.ones((1,number_of_subsets)),np.array([np.mean(temp_minmax_table,axis=1)]).transpose())
                    # create the full table
                    # note: this table has initially the same shape as the observable table
                    temp_full_table = np.column_stack((
                            temp_ids_table,
                            temp_minmax_table,
                            temp_initial_table,
                            temp_estimated_table))
                    # take out unique rows and the inverse vector recreating the rows of the observable table
                    # the inverse vector is critical as it relates the positions in the observable table to positions in each parameter table
                    temp_full_table,temp_position_in_observable_table = np.unique(temp_full_table,axis=0,return_inverse=True)
                    # set the last colum of the full table to nan (estimated value unavailable now)
                    temp_full_table[:,-number_of_subsets:] = np.nan
                    # append the table
                    parameter_tables.append(temp_full_table)
                    # include the inverse vector in the observable data table at the correct position
                    parameter_position_table[:,parameter_idx] = temp_position_in_observable_table
                
                
                line_number_string = ['row']

                ### SAVING OBSERVABLE AND PARAMETER TABLES
                # select formatting for the output tables
                curr_delimiter = '\t'
                curr_precision = 3
                curr_column_width = 20
                all_column_groups = [line_number_string,identifier_names,parameter_position_names,parameter_property_names,observable_names,parameter_names]
                all_data_types = ['d','d','d','f','f','f']
                
                
                # # save parameter and observable tables
                # for parameter_idx,parameter_name in enumerate(parameter_names):
                #     curr_format,curr_header = get_fmt_and_header(parameter_table_columns[parameter_idx],all_column_groups,all_data_types,curr_delimiter,curr_precision,curr_column_width)
                #     np.savetxt(os.path.join(
                #             temp_agb_folder,
                #             'parameter_{}_table_block_{}.txt'.format(parameter_name,current_block_index)),
                #         parameter_tables[parameter_idx],
                #         fmt=curr_format,
                #         delimiter=curr_delimiter,
                #         header=curr_header,
                #         comments='')
                
                curr_format,curr_header = get_fmt_and_header(line_number_string+identifier_names+observable_names+parameter_position_names,all_column_groups,all_data_types,curr_delimiter,curr_precision,curr_column_width)
                np.savetxt(os.path.join(
                        temp_agb_folder,
                        'observable_data_table_block_{}.txt'.format(current_block_index)),
                    np.column_stack((np.arange(identifier_table.shape[0]),identifier_table,observable_table,parameter_position_table)),
                    fmt=curr_format,
                    delimiter=curr_delimiter,
                    header=curr_header,
                    comments='')
            
                
                        # %%
            except Exception as e:
                logging.error('AGB: error during data sampling and tabulating.' + str(e), exc_info=True)
                raise
            
            try:
                # %%
                
                ### CREATE CALIBRATION AND ESTIMATION SUBSETS
                    
                # select rows with available agb information as calibration data and those without as estimation data
                calibration_rows = np.where(np.all(~np.isnan(observable_table),axis=1))[0]
                estimation_rows = np.where(np.any(np.isnan(observable_table),axis=1))[0]
                calibration_sample_ids = np.unique(identifier_table[calibration_rows,0])
                estimation_sample_ids = np.unique(identifier_table[estimation_rows,0])
                
                # calculate subset sizes
                estimation_subset_size = np.int32(np.ceil(len(estimation_sample_ids)/100*proc_conf.AGB.fraction_of_roi_per_test))
                calibration_subset_size = np.int32(np.ceil(len(calibration_sample_ids)/100*proc_conf.AGB.fraction_of_cal_per_test))
                estimation_subset_size = np.int32(np.ceil(len(estimation_sample_ids)/100*80))
                calibration_subset_size = np.int32(np.ceil(len(calibration_sample_ids)/100*80))
                
                # find random data subsetting vectors making sure that the number of calibration and estimation areas
                # is the same in all
                subset_indexing_vectors = []
                number_of_accepted_subsets = 0
                while number_of_accepted_subsets < number_of_subsets:
                    # create a random subset of calibration and estimation samples
                    current_random_estimation_subset = np.sort(np.random.permutation(estimation_sample_ids)[:estimation_subset_size])
                    current_random_calibration_subset = np.sort(np.random.permutation(calibration_sample_ids)[:calibration_subset_size])
                    
                    # calculate the minimal number of calibration and estimation samples for the space-invariant parameters
                    # (for the latter, we use the column with parameter positions in parameter tables - the same value indicates the same parameter)
                    current_calibration_rows = np.isin(identifier_table[:,0],current_random_calibration_subset)
                    current_estimation_rows = np.isin(identifier_table[:,0],current_random_estimation_subset)
                    current_parameter_position_columns = np.where(~np.row_stack(parameter_variabilities)[:,0])[0]
                    min_number_of_calibration_measurements_per_space_invariant_parameter = np.inf
                    min_number_of_estimation_measurements_per_space_invariant_parameter = np.inf
                    # loop through columns with parameter positions
                    for column_idx in current_parameter_position_columns:
                        # calculate the minimal number of samples for all parameter values within this column and all different parameters until the current one
                        min_number_of_calibration_measurements_per_space_invariant_parameter\
                            = np.minimum(min_number_of_calibration_measurements_per_space_invariant_parameter,
                                                            np.min(
                                                                np.unique(
                                                                 parameter_position_table[current_calibration_rows,column_idx],return_counts=True)[1]))
                        min_number_of_estimation_measurements_per_space_invariant_parameter\
                            = np.minimum(min_number_of_estimation_measurements_per_space_invariant_parameter,
                                                            np.min(
                                                                np.unique(
                                                                 parameter_position_table[current_estimation_rows,column_idx],return_counts=True)[1]))
                    # if the minimal number of samples is larger than the one specified in the xml configuration file, accept this subset
                    # (at the moment, we don't perform other tests, which means that subsets may be repeated)
                    if (min_number_of_calibration_measurements_per_space_invariant_parameter>proc_conf.AGB.min_number_of_cals_per_test) & \
                        (min_number_of_estimation_measurements_per_space_invariant_parameter>proc_conf.AGB.min_number_of_rois_per_test):
                            subset_indexing_vectors.append(np.isin(identifier_table[:,0],np.sort(np.concatenate((current_random_calibration_subset,current_random_estimation_subset)))))
                            number_of_accepted_subsets += 1
                            
                
                
                
                
                
                ### RUN INVERSION FOR THE CURRENT BLOCK AND TABLE USING ALL SELECTED SUBSETS
            
            
            
                # loop through calibration subsets
                number_of_successful_tests = 0
                for current_subset in subset_indexing_vectors:
                    


                                            
                    # take out a subtable with indices for each parameter in the output parameter tables
                    current_parameter_position_table = parameter_position_table[current_subset,:]
                    # take out a subtable with observables
                    current_observable_table = observable_table[current_subset,:]
                    
                        
                    
                    # convert these positions to positions in a single vector x, where unique
                    # parameter values are ordered first by rows (samples) and then by columns (independent parameters)
                    # this requires renumbering and introducing offsets for each column
                    offset = 0
                    current_parameter_index_in_beta_table = [] # position in a single x-vector
                    current_parameter_position_lut = [] # lut for converting between parameter id and column and parameter position in x
                    for parameter_idx in range(number_of_parameters):
                        # add -1 at the beginning to avoid vectors with one element (which will not work with interp1d)
                        old_indices = np.concatenate((-1*np.ones(1),np.unique(current_parameter_position_table[:,parameter_idx])))
                        # new indices is a simple sequence from 0 to number of parameters-1 + offset due to previous parameters
                        new_indices = np.arange(len(old_indices))-1+offset
                        # convert parameter indices and add to the list
                        current_parameter_index_in_beta_table.append(sp.interpolate.interp1d(old_indices,new_indices,kind='nearest')(current_parameter_position_table[:,parameter_idx]))
                        # save the lut, removing the first, unnecessary element
                        current_parameter_position_lut.append(
                            np.column_stack((new_indices[old_indices>-1],
                                             parameter_idx+np.zeros(len(old_indices[old_indices>-1])),
                                             old_indices[old_indices>-1],
                                             )))
                        # update offset based on current parameter column
                        offset = np.max(new_indices)+1
                    # convert the list of vectors to an array
                    current_parameter_index_in_beta_table = np.int32(np.column_stack(current_parameter_index_in_beta_table))
                    # stack all luts to one lut
                    current_parameter_position_lut = np.row_stack(current_parameter_position_lut)
                    # length of beta vector
                    length_of_beta_vector = current_parameter_position_lut.shape[0]
                    
                    # create a table of min, max, and initial parameter values (requires looping through parameter tables and extracting relevant rows)
                    beta_min_max_initial_table = np.nan * np.zeros((length_of_beta_vector,3))
                    for parameter_idx in np.int32(np.unique(current_parameter_position_lut[:,1])):
                        current_rows_in_lut = current_parameter_position_lut[:,1]==parameter_idx
                        current_positions_in_parameter_table = np.int32(current_parameter_position_lut[current_rows_in_lut,2])
                        beta_min_max_initial_table[current_rows_in_lut,:] = \
                            parameter_tables[parameter_idx][current_positions_in_parameter_table,-(number_of_subsets+3):-number_of_subsets]
                    # extract min, max, intiial values
                    beta_low = beta_min_max_initial_table[:,0]
                    beta_high = beta_min_max_initial_table[:,1]
                    beta_initial = beta_min_max_initial_table[:,2]
                    # x_initial = parameter_transfer_function(beta_initial,beta_low,beta_high,True)
                    x_initial = np.pi/4*np.random.rand(length_of_beta_vector)
                    
                    cal_rows = np.all(~np.isnan(current_observable_table),axis=1)
                    est_rows = np.any(np.isnan(current_observable_table),axis=1)
                    
                    
                    calibration_formula = convert_formula(formula,'calibration_data_table',observable_names + parameter_names,0)
                    estimation_formula = convert_formula(formula,'estimation_data_table',observable_names + parameter_names,1)
                        
                    # calibration_data_table = np.column_stack((current_observable_table[cal_rows,:],beta_initial[current_parameter_index_in_beta_table[cal_rows,:]]))
                    
                    def J(x,args):
                        beta_low,beta_high,parameter_transfer_function,current_observable_table,cal_rows,est_rows,current_parameter_index_in_beta_table,calibration_formula,estimation_formula = args
                        beta = parameter_transfer_function(x,beta_low,beta_high,False)
                        
                        calibration_data_table = np.column_stack((current_observable_table[cal_rows,:],beta[current_parameter_index_in_beta_table[cal_rows,:]]))
                        estimation_data_table = np.column_stack((current_observable_table[est_rows,:],beta[current_parameter_index_in_beta_table[est_rows,:]]))
                        # J = np.mean(eval(estimation_formula))+
                        J = np.sqrt(np.mean(eval(calibration_formula))/2+np.mean(eval(estimation_formula))/2)
                        # print(J)
                        return J
                    # j = J(x_initial,args)
                    
                    args = [beta_low,beta_high,parameter_transfer_function,current_observable_table,cal_rows,est_rows,current_parameter_index_in_beta_table,calibration_formula,estimation_formula]
                    mod = sp.optimize.minimize(J,x_initial,args,method='BFGS')
                    
                    x_estimated = mod.x
                    x_estimated[np.abs(x_initial-x_estimated)<1e-6] = np.nan
                    beta = parameter_transfer_function(x_estimated, beta_low,beta_high,False)
                    print(mod.message,J(x_estimated,args))
                    if mod.success:
                        number_of_successful_tests += 1
                        for parameter_idx in range(number_of_parameters):
                            curr_rows = current_parameter_position_lut[:,1]==parameter_idx
                            
                            # error
                            parameter_tables[parameter_idx][np.int32(current_parameter_position_lut[curr_rows,2]),-number_of_subsets+number_of_successful_tests-1] = \
                                beta[np.int32(current_parameter_position_lut[curr_rows,0])]
                            # parameter_table[:,parameter_idx] = np.nanmean(parameter_tables[parameter_idx][:,-number_of_subsets:],axis=1)[np.int32(parameter_position_table[:,parameter_idx])]
                
                parameter_table = np.nan * np.zeros((number_of_rows_in_observable_table, number_of_parameters))    
                for parameter_idx,parameter_name in enumerate(parameter_names):
                # for parameter_idx in range(number_of_parameters):
                    curr_rows = current_parameter_position_lut[:,1]==parameter_idx
                    
                    # error
                    # parameter_tables[parameter_idx][np.int32(current_parameter_position_lut[curr_rows,2]),-number_of_subsets+number_of_successful_tests-1] = \
                    #     beta[np.int32(current_parameter_position_lut[curr_rows,0])]
                    parameter_table[:,parameter_idx] = np.nanmean(parameter_tables[parameter_idx][:,-number_of_subsets:],axis=1)[np.int32(parameter_position_table[:,parameter_idx])]
            
                    curr_format,curr_header = get_fmt_and_header(np.concatenate((np.array(line_number_string),parameter_table_columns[parameter_idx])),all_column_groups,all_data_types,curr_delimiter,curr_precision,curr_column_width)
                    np.savetxt(os.path.join(
                            temp_agb_folder,
                            'parameter_{}_table_block_{}.txt'.format(parameter_name,current_block_index)),
                        np.column_stack((np.arange(parameter_tables[parameter_idx].shape[0]),parameter_tables[parameter_idx])),
                        fmt=curr_format,
                        delimiter=curr_delimiter,
                        header=curr_header,
                        comments='')
                
                curr_format,curr_header = get_fmt_and_header(line_number_string+identifier_names+observable_names+parameter_names,all_column_groups,all_data_types,curr_delimiter,curr_precision,curr_column_width)
                np.savetxt(os.path.join(
                        temp_agb_folder,
                        'final_data_table_block_{}.txt'.format(current_block_index)),
                    np.column_stack((np.arange(identifier_table.shape[0]),identifier_table,observable_table,parameter_table)),
                    fmt=curr_format,
                    delimiter=curr_delimiter,
                    header=curr_header,
                    comments='')
                            
                    # # %%
                    
                    # columns_with_agb_parameters = np.where(np.array([curr_source.split('_')[0]=='agb' for curr_source in parameter_names]))[0]
                    # columns_with_sar_parameters = np.where(np.array([curr_source.split('_')[0]!='agb' for curr_source in parameter_names]))[0]
                    
                    # # convert these positions to positions in a single vector x, where unique
                    # # parameter values are ordered first by rows (samples) and then by columns (independent parameters)
                    # # this requires renumbering and introducing offsets for each column
                    # offset = 0
                    # current_parameter_index_in_beta_table = [] # position in a single x-vector
                    # current_parameter_position_lut = [] # lut for converting between parameter id and column and parameter position in x
                    # for parameter_idx in columns_with_sar_parameters:
                    #     # add -1 at the beginning to avoid vectors with one element (which will not work with interp1d)
                    #     old_indices = np.concatenate((-1*np.ones(1),np.unique(current_parameter_position_table[:,parameter_idx])))
                    #     # new indices is a simple sequence from 0 to number of parameters-1 + offset due to previous parameters
                    #     new_indices = np.arange(len(old_indices))-1+offset
                    #     # convert parameter indices and add to the list
                    #     current_parameter_index_in_beta_table.append(sp.interpolate.interp1d(old_indices,new_indices,kind='nearest')(current_parameter_position_table[:,parameter_idx]))
                    #     # save the lut, removing the first, unnecessary element
                    #     current_parameter_position_lut.append(
                    #         np.column_stack((new_indices[old_indices>-1],
                    #                          parameter_idx+np.zeros(len(old_indices[old_indices>-1])),
                    #                          old_indices[old_indices>-1],
                    #                          )))
                    #     # update offset based on current parameter column
                    #     offset = np.max(new_indices)+1
                    # # convert the list of vectors to an array
                    # current_parameter_index_in_beta_table = np.int32(np.column_stack(current_parameter_index_in_beta_table))
                    # # stack all luts to one lut
                    # current_parameter_position_lut = np.row_stack(current_parameter_position_lut)
                    # # length of beta vector
                    # length_of_beta_vector = current_parameter_position_lut.shape[0]
                    
                    # # create a table of min, max, and initial parameter values (requires looping through parameter tables and extracting relevant rows)
                    # beta_min_max_initial_table = np.nan * np.zeros((length_of_beta_vector,3))
                    # for parameter_idx in np.int32(np.unique(current_parameter_position_lut[:,1])):
                    #     current_rows_in_lut = current_parameter_position_lut[:,1]==parameter_idx
                    #     current_positions_in_parameter_table = np.int32(current_parameter_position_lut[current_rows_in_lut,2])
                    #     beta_min_max_initial_table[current_rows_in_lut,:] = \
                    #         parameter_tables[parameter_idx][current_positions_in_parameter_table,-(number_of_subsets+3):-number_of_subsets]
                    # # extract min, max, intiial values
                    # beta_low = beta_min_max_initial_table[:,0]
                    # beta_high = beta_min_max_initial_table[:,1]
                    # beta_initial = beta_min_max_initial_table[:,2]
                    # # x_initial = parameter_transfer_function(beta_initial,beta_low,beta_high,True)
                    # x_initial = np.pi/4*np.random.rand(length_of_beta_vector)
                    
                    
                    
                    # new_formula = convert_formula(formula,'new_data_table',observable_names + parameter_names,1)
                   
                    
                    # # calibration_data_table = np.column_stack((current_observable_table[cal_rows,:],beta_initial[current_parameter_index_in_beta_table[cal_rows,:]]))
                    
                    # def J(x,args):
                    #     beta_low,beta_high,parameter_transfer_function,current_observable_table,valid_rows,current_parameter_index_in_beta_table,new_formula = args
                    #     beta = parameter_transfer_function(x,beta_low,beta_high,False)
                        
                    #     new_data_table = np.column_stack((current_observable_table[valid_rows,:],beta[current_parameter_index_in_beta_table[valid_rows,:]]))
                        
                    #     # J = np.mean(eval(estimation_formula))+
                    #     J = np.sqrt(np.mean(eval(new_formula)))
                    #     # print(J)
                    #     return J
                    
                    # # j = J(x_initial,args)
                    
                    # args = [beta_low,beta_high,parameter_transfer_function,np.column_stack((current_observable_table,,cal_rows,est_rows,current_parameter_index_in_beta_table,calibration_formula,estimation_formula]
                    # mod = sp.optimize.minimize(J,x_initial,args,method='BFGS')
                    
                    # x_estimated = mod.x
                    # x_estimated[np.abs(x_initial-x_estimated)<1e-6] = np.nan
                    # beta = parameter_transfer_function(x_estimated, beta_low,beta_high,False)
                    # print(mod.message,J(x_estimated,args))
                    # # %%
                
                
                
                ### now we average the estimated AGBs
                
                
                
                # save parameter and observable tables
                
                
                # %%
                    # estimation_formula = convert_formula(formula,'data_table',observable_names + parameter_names,1)
                    
                    
                    
                    
                    
                    # def estimate_for_table(x_initial):
                        
                    #     def J(x):
                                
                    #         # take out the rows with calibration data
                    #         beta = parameter_transfer_function(x,beta_low,beta_high,False)
                    #         data_table = np.column_stack((current_observable_table[cal_rows,:],beta[current_parameter_index_in_beta_table[cal_rows,:]]))
                    #         cal_sum = np.mean(eval(calibration_formula))
                    #         # est_sum = np.mean(eval(estimation_formula))
                    #         data_table = np.column_stack((current_observable_table[est_rows,:],beta[current_parameter_index_in_beta_table[est_rows,:]]))
                    #         # cal_sum = np.mean(eval(calibration_formula))
                    #         est_sum = np.mean(eval(estimation_formula))
                    #         return cal_sum+est_sum
                    #     return sp.optimize.minimize(J,x_initial)
                    
                    # estimate_for_table(x_initial)
                    
                    
                    # error
                    # # print(x0)
                    # # error
                    # def estimate(x_initial,
                    #                                 beta_low,
                    #                                 beta_high,
                    #                                 parameter_transfer_method,
                    #                                 current_estimation_parameter_position_in_beta_vector,
                    #                                 current_calibration_parameter_position_in_beta_vector,
                    #                                 current_estimation_observable_table,
                    #                                 current_calibration_observable_table,
                    #                                 estimation_formula_terms,
                    #                                 calibration_formula_terms):
                    #     def J(x,beta_low,beta_high,parameter_transfer_method,current_parameter_index_in_beta_table,current_observable_table,formula_terms):
                            
                    #         def parameter_transfer_function(x,pmin,pmax,kind,x_is_p=False):
                    #             # note: no check of x, xmax, xmin and kind is done here, 
                    #             # it is assumed that the inputs are correct
                    #             if kind=='sin2' :
                    #                 if not x_is_p:
                    #                     return pmin+(pmax-pmin)*np.sin(x)**2
                    #                 else:
                    #                     return np.arcsin(np.sqrt((x-pmin)/(pmax-pmin)))
                    #         def apply_formula(table,formula):
                    #             total = 0
                    #             for independent_term in formula:
                    #                 curr_sum = 0
                    #                 for added_term in independent_term:
                    #                     curr_product = 1
                    #                     for multiplied_term in added_term:
                    #                         curr_product *= table[:,multiplied_term] 
                    #                         # print(table[-15,multiplied_term] )
                    #                     curr_sum += curr_product
                    #                 total += (curr_sum)**2
                    #             return total
                    #         beta = parameter_transfer_function(x,beta_low,beta_high,parameter_transfer_method)
                    #         current_parameter_table = beta[current_parameter_index_in_beta_table]
                    #         current_combined_table = np.column_stack((current_observable_table,current_parameter_table))
                            
                    #         # print(np.any(np.isnan(calibration_table),axis=0))
                    #         # return apply_formula(estimation_table,parsed_formula_estimation)#+
                    #         # print(current_combined_table)
                    #         # tot = np.nanmean(apply_formula(current_combined_table,estimation_formula_terms))+\
                    #         tot = np.mean(apply_formula(current_combined_table,formula_terms))
                    #         print(tot)
                    #         return tot
                                
                            
                        
                        
                    #     return sp.optimize.minimize(lambda x:J(x,
                    #                                     beta_low,
                    #                                     beta_high,
                    #                                     parameter_transfer_method,
                    #                                     current_calibration_parameter_position_in_beta_vector,
                    #                                     current_calibration_observable_table,
                    #                                     calibration_formula_terms),x_initial,method='Nelder-Mead')
                        
                            
                    # print(estimate(x_initial,
                    #                             beta_low,
                    #                             beta_high,
                    #                             parameter_transfer_method,
                    #                             current_parameter_index_in_beta_table[current_estimation_rows,:],
                    #                             current_observable_table,
                    #                             estimation_formula_terms,
                    #                             calibration_formula_terms))
                    # error
                    #     # print(current_parameter_index_in_beta_table,current_parameter_position_table)
                        
                    #     # print(np.row_stack([np.array([current_position_in_x,np.where([current_parameter_index_in_beta_table==current_position_in_x])[0][0],
                    #     #           np.where([current_parameter_index_in_beta_table==current_position_in_x])[1][0]]) for current_position_in_x in np.unique(current_parameter_index_in_beta_table)]))
                        
                    #     # for parameter_position_column,parameter_table in zip(current_parameter_position_table.transpose(),parameter_tables):
                    #     #     print(parameter_table[np.int32(parameter_position_column),-4])
                        
                    #     # print(current_parameter_index_in_beta_table)
                    #     # # merge parameter tables
                    #     # parameter_tables_merged = np.row_stack([curr_parameter_table[:,-4:-1] for curr_parameter_table in parameter_tables])
                    #     # # calculate index offsets for first element of each table
                    #     # parameter_index_offsets = [curr_parameter_table.shape]
                    #     # p_low = p_rands[:,0]
                    #     # p_high = p_rands[:,1]
                    #     # p_init = p_rands[:,2]
                        
        
                    #     # print(parameter_transfer_function(p_init,p_low,p_high,parameter_transfer_method,True))
                    #     # x_init = 
                        
                            
                        
                    
                    # # return []
                
                
                

                # %%
            
            
            # all_est_sample_ids = np.unique(observable_est_data_table[:,0])
            # all_cal_sample_ids = np.unique(observable_cal_data_table[:,0])
            
            
            # # if (est_subset_size<proc_conf.AGB.)
            
            # # creating subsets of calibration and estimation data
            # number_of_accepted_tests = 0
            # est_subset_ids = []
            # cal_subset_ids = []
            # while number_of_accepted_tests < number_of_subsets:
            #     temp_est_subset = all_est_sample_ids[np.random.permutation(np.arange(len(all_est_sample_ids)))[:est_subset_size]]
            #     temp_cal_subset = all_cal_sample_ids[np.random.permutation(np.arange(len(all_cal_sample_ids)))[:cal_subset_size]]
            #     # est_subset_ids.append()
            #     # cal_subset_ids.append()
                
            #     temp_parameter_indices_in_cal_data = observable_cal_data_table[np.isin(observable_cal_data_table[:,0],temp_cal_subset),-number_of_parameters:]
                
            #     for parameter_idx in range(number_of_parameters):
            #         if not parameter_variabilities[parameter_idx][0]:
            #             print(np.min(np.unique(temp_parameter_indices_in_cal_data[:,parameter_idx],return_counts=True)[1]))
            #     number_of_accepted_tests += 1
            # np.unique(observable_cal_data_table[:,-number_of_parameters],return_counts=True)
            
            
# 0 1           # %%
            except Exception as e:
                logging.error('AGB: error during parameter estimation from tables.' + str(e), exc_info=True)
                raise
            
            
            
            
            # %%
            # temp_parameter_positions = np.column_stack(position_ids)
            
            
            
            
            
            
            
            
            
        
            
            
               
            # np.save(
            #     os.path.join(
            #         temp_agb_folder,
            #         'observable_table_block_{}'.format(current_block_index),
            #     ),
            #     observable_table,
            # )
            #     # temp_agb_vec = np.nan * np.zeros( (len(temp_sample_id),number_of_observables))
            #     # # rearrange agb data into tables
            #     # for descriptor_idx,agb_descriptor in enumerate(agb_descriptors):
            #     #     temp_data_vec[:,observable_idx] = observable_data_sampled[:,observable_idx,:].flatten()
                
            # # rearrange agb descriptors into tables
            # for descriptor_idx,descriptor_name in enumerate(agb_descriptors):
            #     np.column_stack(
            #         (
            #             sample_id_array.flatten(),
            #             forest_class_sampled.flatten(),
            #             agb_data_sampled[:,descriptor_idx]))
                

                            # # # Insert in table
                            # # sigma_tab[counter_stack_pols, :] = data_roi_means_vec

                            # progressive_stack_idx_vec = (
                            #     np.ones(number_of_samples_on_grid, dtype='int') * progressive_stack_idx
                            # )
                            # progressive_stack_idx_vec[
                            #     np.isnan(data_roi_means_vec)
                            # ] = -1  # invalid value

                            # acqid_data_sampled[0,, :] = progressive_stack_idx

                            # progressive_stack_idx_vec = (
                            #     np.ones(number_of_pixels, dtype='int') * progressive_stack_idx
                            # )
                            # progressive_stack_idx_vec[
                            #     np.isnan(stack_data_interp.flatten())
                            # ] = -1  # invalid value

                            # acqid_tab_pixels[counter_stack_pols, :] = progressive_stack_idx_vec

                            # counter_stack_pols = counter_stack_pols + 1

                        

                        # # cycle over different angles to read
                        # for angle_idx,angle_name in enumerate(angle_names):
                        #     # cycle over all the equi7 tiles, interpolate over pixel grid and mean them togheter
                            
                        #     angle_eq7_parent_dir_name = os.path.join(
                        #         stack_path, angle_name, equi7_grid_name
                        #     )
    
                        #     stack_angle_interp = np.NaN * np.zeros(
                        #         (len(pixel_axis_north), len(pixel_axis_east)), dtype='float'
                        #     )
                        #     for equi7_tile_name in equi7_tiles_names:
    
                        #         angle_tiff_name = os.path.join(
                        #             angle_eq7_parent_dir_name,
                        #             equi7_tile_name, 
                        #             angle_name
                        #             + '_'
                        #             + os.path.basename(stack_path)
                        #             + '_'
                        #             + equi7_grid_name[6:]
                        #             + '_'
                        #             + equi7_tile_name
                        #             + '.tif',
                        #         )
    
                        #         stack_angle_interp_curr = interp2d_wrapper(
                        #             angle_tiff_name, 1, pixel_axis_east, pixel_axis_north, fill_value=np.NaN
                        #         )
    
                        #         # mean all the equi7 tiles
                        #         stack_angle_interp = merge_agb_intermediate(
                        #             stack_angle_interp, stack_angle_interp_curr, method='nan_mean'
                        #         )
    
                        #         projection_curr = get_projection_from_path(angle_tiff_name)
                        #         if projection_prev and (projection_prev != projection_curr):
                        #             err_str = 'multiple equi7 projections are not supported.'
                        #             logging.error(err_str)
                        #             raise ImportError(err_str)
                        #         projection_prev = projection_curr

                        #     # masking the stack:
                        #     stack_angle_interp[forest_class_map_interp == 0] = np.NaN
    
                        #     # theta_tab_pixels[counter_stacks, :] = stack_theta_interp.flatten()
    
                        #     # Mean on ROI
                        #     angle_data_sampled[:,angle_idx,stack_idx] = mean_on_rois(
                        #         stack_angle_interp,
                        #         pixel_mesh_east,
                        #         pixel_mesh_north,
                        #         sampling_axis_east,
                        #         sample_spacing_east,
                        #         sampling_axis_north,
                        #         sample_spacing_north,
                        #         angle_averaging_methods[angle_idx],
                        #     )
    
                            # # Insert in table
                            # theta_tab[counter_stacks, :] = theta_roi_means_vec
    
                            # counter_stacks = counter_stacks + 1

                # # replicate equal values for all the 3 polarizations:
                # theta_tab = np.kron(theta_tab, np.ones((3, 1)))
                # theta_tab_pixels = np.kron(theta_tab_pixels, np.ones((3, 1)))

                # N_laypol = sigma_tab.shape[0]  # number of stacks x number of polarizations'
                # N_lay = np.int32(N_laypol / 3)  # number of stacks
                # create polarisation ID table matching sigma data table in shape and size
                # polid_tab = np.kron(
                #     np.ones((N_lay, 1)),
                #     np.kron(np.array([np.arange(3)]).transpose(), np.ones((1, number_of_samples_on_grid))),
                # )
                # polid_tab_pixels = np.kron(
                #     np.ones((N_lay, 1)),
                #     np.kron(np.array([np.arange(3)]).transpose(), np.ones((1, number_of_pixels))),
                # )

                # # save and delete all the pixel tabs:
                # np.save(
                #     os.path.join(
                #         sigma_tab_pixel_folder,
                #         'sigma_tab_pixels_block_{}'.format(current_block_index),
                #     ),
                #     sigma_tab_pixels,
                # )
                # np.save(
                #     os.path.join(
                #         acq_tab_pixel_folder, 'acq_tab_pixels_block_{}'.format(current_block_index)
                #     ),
                #     acqid_tab_pixels,
                # )
                # np.save(
                #     os.path.join(
                #         theta_tab_pixel_folder,
                #         'theta_tab_pixels_block_{}'.format(current_block_index),
                #     ),
                #     theta_tab_pixels,
                # )
                # np.save(
                #     os.path.join(
                #         polid_tab_pixel_folder,
                #         'polid_tab_pixels_block_{}'.format(current_block_index),
                #     ),
                #     polid_tab_pixels,
                # )
                # del sigma_tab_pixels, acqid_tab_pixels, theta_tab_pixels, polid_tab_pixels

                


                # ## FILTER OUT INVALID ACQUISITIONS
                
                # # set all measurements to nan for samples that do not fulfill certain requirements on sigma, theta, etc.
                # for stack_idx in range(number_of_stacks):
                #     for obs_idx,observable_name in enumerate(observable_names):
                #         if observable_name.split('_')[0] == 'sigma':
                #             invalid_rows = ( 
                #                     observable_data_sampled[:,obs_idx,stack_idx]<=0
                #                     ) | (
                #                         ~np.isfinite(observable_data_sampled[:,obs_idx,stack_idx])
                #                     ) | (
                #                         np.isnan(observable_data_sampled[:,obs_idx,stack_idx])
                #                     )
                #             observable_data_sampled[invalid_rows,:,stack_idx] = np.nan
                #         elif observable_name.split('_')[0] == 'theta':
                #             invalid_rows = ( 
                #                     observable_data_sampled[:,obs_idx,stack_idx]<0
                #                     ) | (
                #                         ~np.isfinite(observable_data_sampled[:,obs_idx,stack_idx])
                #                     ) | (
                #                         np.isnan(observable_data_sampled[:,obs_idx,stack_idx])
                #                     )
                #             observable_data_sampled[invalid_rows,:,stack_idx] = np.nan
                            
                            
                # unique_acqs = acqid_data_sampled.flatten()
                # acq_count = np.sum(np.all(~np.isnan(observable_data_sampled),axis=1),axis=0)
                # # find independent measurements within current tile and count them
                # # unique_acqs, acq_count = np.unique(acqid_tab, return_counts=True)
                # # extract acquisitions meeting the requirements on minimal number of independent measurements
                # valid_acqs = np.isin(
                #     stack_info_table[:, 0],
                #     unique_acqs[acq_count >= proc_conf.AGB.min_number_of_rois_per_stack],
                # )
                # # mask out invalid acquisitions (either too few measurements or there are nans in sigma or theta images)
                # # acqid_tab[~np.isin(acquisition_, stack_info_table[valid_acqs, 0]) ] = -1

                # # # mask out data in sigma0 and local tables
                # # sigma_tab[acqid_tab == -1] = np.nan
                # # theta_tab[acqid_tab == -1] = np.nan


            # # make sure we have enough data
            # if not sum(valid_acqs == True):
            #     logging.info(
            #         'skipping block #{} due to too invalid data points'.format(current_block_index)
            #     )
            #     # swap the flag
            #     block_finished_flag[current_block_index] = True
            #     # remove current par block from list
            #     block_order = block_order[block_order != current_block_index]
            #     # select next par block
            #     if len(block_order) > 0:
            #         # if there is at least one left, just take the next closest to CALdata
            #         current_block_index = block_order[0]
            #         continue
            #     else:
            #         break

            # ## PREPARE PARAMETER ID TABLES
            # # create roi id table
            # roiid_tab = np.int32(np.array([np.arange(number_of_samples_on_grid)]) * np.ones((N_laypol, 1)))
            # # create look-up tables for acquisition ID and unique parameter ID
            # par_id_luts = [
            #     np.unique(
            #         stack_info_table[valid_acqs, :] * np.array([curr_var]), axis=0, return_inverse=True
            #     )[1]
            #     for curr_var in parameter_variabilities
            # ]
            # # create parameter id tables matching data tables
            # # stack_info_table is ordered differently from other tables with respect to acquisition id and polarization
            # # tableLookupInt function takes that into account
            # par_tabs = [
            #     tableLookupInt(
            #         stack_info_table[valid_acqs, :2],
            #         curr_lut,
            #         np.column_stack((polid_tab.flatten(), acqid_tab.flatten())),
            #     ).reshape(acqid_tab.shape)
            #     for curr_lut in par_id_luts
            # ]
            # # extract the number of unique values for each parameter
            # all_N_par = np.int32(
            #     np.array([np.nanmax(par_tab) + 1 for par_tab in par_tabs + [roiid_tab]])
            # )
            # # calculate the offset from first parameter in the large vector
            # par_offsets = np.int32(np.concatenate((np.zeros(1), np.cumsum(all_N_par))))[:-1]

            # if at least one tile has been finished, load the low-resolution estimates for the relevant pixels
            
            # # replicate N_laypol times
            # cal_tab = np.kron(cal_data_roi_means_vec, np.ones((N_laypol, 1)))
            # std1_tab_db = np.kron(cal_std1_roi_means_vec, np.ones((N_laypol, 1)))
            # std2_tab_db = np.kron(cal_std2_roi_means_vec, np.ones((N_laypol, 1)))
            # std3_tab_db = np.kron(cal_std3_roi_means_vec, np.ones((N_laypol, 1)))

            ## PREPARE CALIBRATION SETS
            # now, lets take out calibration ROIs and randomise two of them to be used for calibration in multiple tests
            all_cal_ids = np.unique(roiid_tab[~np.isnan(cal_tab) & (acqid_tab > -1)])
            all_roi_ids = np.setdiff1d(
                np.unique(roiid_tab[acqid_tab > -1]), all_cal_ids
            )  # tolgo le cal dalle roi
            # number of CALs in a subset (approximately half the total number of CALs)
            N_cal_sub = np.int32(
                np.floor(proc_conf.AGB.fraction_of_cal_per_test / 100 * len(all_cal_ids))
            )
            N_roi_sub = np.int32(
                np.floor(proc_conf.AGB.fraction_of_roi_per_test / 100 * len(all_roi_ids))
            )

            logging.info('Number of cals in current block :{}'.format(N_cal_sub))
            logging.info('Number of rois in current block :{}'.format(N_roi_sub))

            # make sure we have enough data
            if not (
                (N_cal_sub >= proc_conf.AGB.min_number_of_cals_per_test)
                & (N_roi_sub >= proc_conf.AGB.min_number_of_rois_per_test)
            ):
                logging.info(
                    'skipping block #{} due to too few data points'.format(current_block_index)
                )
                # swap the flag
                block_finished_flag[current_block_index] = True
                # remove current par block from list
                block_order = block_order[block_order != current_block_index]
                # select next par block
                if len(block_order) > 0:
                    # if there is at least one left, just take the next closest to CALdata
                    current_block_index = block_order[0]
                    continue
                else:
                    break

            # create random subsets of N_cal and N_roi
            cal_sets = all_cal_ids[
                np.unique(
                    np.row_stack(
                        [
                            np.sort(np.random.permutation(np.arange(len(all_cal_ids)))[:N_cal_sub])
                            for ii in np.arange(10000)
                        ]
                    ),
                    axis=0,
                )
            ]
            cal_sets = cal_sets[np.random.permutation(np.arange(cal_sets.shape[0]))[:number_of_subsets], :]
            roi_sets = all_roi_ids[
                np.unique(
                    np.row_stack(
                        [
                            np.sort(np.random.permutation(np.arange(len(all_roi_ids)))[:N_roi_sub])
                            for ii in np.arange(10000)
                        ]
                    ),
                    axis=0,
                )
            ]
            roi_sets = roi_sets[np.random.permutation(np.arange(roi_sets.shape[0]))[:number_of_subsets], :]
            # allocate parameter table
            curr_output_table = np.nan * np.zeros((N_laypol, number_of_samples_on_grid, 6, number_of_subsets))
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
                eps1 = np.random.randn(1, number_of_samples_on_grid)
                eps2 = np.random.randn(N_laypol, 1)
                eps3 = np.random.randn()
                # extract data for rois and cals
                s_cal = forward(sigma_tab[curr_calmask])
                c_cal = forward(np.cos(theta_tab[curr_calmask]))
                if proc_conf.AGB.add_variability_on_cal_data:
                    w_cal = (
                        forward(np.maximum(1, cal_tab[curr_calmask]))
                        + (std1_tab_db * eps1)[curr_calmask]
                        + (std2_tab_db * eps2)[curr_calmask]
                        + (std3_tab_db * eps3)[curr_calmask]
                    )
                else:
                    w_cal = forward(np.maximum(1, cal_tab[curr_calmask]))

                s_roi = forward(sigma_tab[curr_roimask])
                c_roi = forward(np.cos(theta_tab[curr_roimask]))
                # number of CAL and ROI measurements
                N = len(s_roi)
                M = len(s_cal)
                # create index vectors for rois and cals
                i0_l_cal, i0_a_cal, i0_n_cal = [
                    par_tab[curr_calmask] + par_offset
                    for par_tab, par_offset in zip(par_tabs, par_offsets)
                ]
                i0_l_roi, i0_a_roi, i0_n_roi = [
                    par_tab[curr_roimask] + par_offset
                    for par_tab, par_offset in zip(par_tabs, par_offsets)
                ]
                i0_w_roi = roiid_tab[curr_roimask] + par_offsets[-1]
                # regularize the indices to get minimal size x-vector
                i_l_cal, i_l_roi, lut_l = regularizeIndices(i0_l_cal, i0_l_roi, 0)
                i_a_cal, i_a_roi, lut_a = regularizeIndices(
                    i0_a_cal, i0_a_roi, np.max(lut_l[:, 1]) + 1
                )
                i_n_cal, i_n_roi, lut_n = regularizeIndices(
                    i0_n_cal, i0_n_roi, np.max(lut_a[:, 1]) + 1
                )
                i_w_roi, i_w_roi, lut_w = regularizeIndices(
                    i0_w_roi, i0_w_roi, np.max(lut_n[:, 1]) + 1
                )
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
                logging.info(
                    '\tMinimisation for test %d/%d finished (%s)'
                    % (i_test + 1, number_of_subsets, mod.message)
                )
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
                    s_roi
                    - (n0 + dn * np.sin(x[i_n_roi]) ** 2) * c_roi
                    - (l0 + dl * np.sin(x[i_l_roi]) ** 2)
                ) * (a0 + da * np.sin(x[i_a_roi]) ** 2)
                curr_output_table[curr_calmask, 4, i_test] = (
                    s_cal
                    - (n0 + dn * np.sin(x[i_n_cal]) ** 2) * c_cal
                    - (l0 + dl * np.sin(x[i_l_cal]) ** 2)
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
                        / np.nansum(
                            curr_output_table[ilay * 3 + np.arange(3), :, 1, :] ** 2, axis=0
                        )
                        for ilay in np.int32(np.arange(N_lay))
                    ]
                )
            )
            # calculate the first residual statistics from the difference between these biomass estimates and the true biomass values for CALs
            res_1 = forward(agbest_tab_eq) - forward(
                np.array([[cal_tab[0, :]]]) * np.ones((number_of_subsets, N_lay, 1))
            ).transpose([1, 2, 0])
            agbstd_1 = np.nanstd(res_1)
            bias_1 = np.nanmean(res_1)
            # re-calculate biomass estimates for each individual triplet and test using the equation and bias estimates
            agbest_tab_eq = inverse(
                (
                    np.array(
                        [
                            np.nansum(curr_output_table[ilay * 3 + np.arange(3), :, 4, :], axis=0)
                            / np.nansum(
                                curr_output_table[ilay * 3 + np.arange(3), :, 1, :] ** 2, axis=0
                            )
                            for ilay in np.int32(np.arange(N_lay))
                        ]
                    )
                    - bias_1
                )
            )
            # calculate biomass estimates averaged over triplets for each individual test using the equation
            agbest_tab_eqbar = inverse(
                (
                    np.nansum(curr_output_table[:, :, 4, :], axis=0)
                    / np.nansum(curr_output_table[:, :, 1, :] ** 2, axis=0)
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
                    np.nansum(curr_output_table[:, :, 4, :], axis=0)
                    / np.nansum(curr_output_table[:, :, 1, :] ** 2, axis=0)
                    - bias_1
                    - np.array([bias_2]).transpose() * np.ones((1, number_of_subsets))
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
            res_3 = forward([agbest]).transpose() * np.ones((1, number_of_subsets)) - forward(
                agbest_tab_eqbar
            )
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
            i_l, i_a, i_n = [
                par_tab[curr_mask] + par_offset
                for par_tab, par_offset in zip(par_tabs, par_offsets)
            ]
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
                            current_block_index,
                            block_corner_coordinates_east[current_block_index],
                            block_corner_coordinates_north[current_block_index],
                            block_corner_coordinates_east[current_block_index] + block_size_east,
                            block_corner_coordinates_north[current_block_index] + block_size_north,
                        ]
                    )
                    * np.ones((np.sum(valid_acqs), 1)),
                    stack_info_table[valid_acqs, 0:2],
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
                agbest.reshape(len(sampling_axis_north), len(sampling_axis_east)),
                (agbest * 0 + agbstd_1).reshape(len(sampling_axis_north), len(sampling_axis_east)),
                (agbest * 0 + agbstd_2).reshape(len(sampling_axis_north), len(sampling_axis_east)),
                (agbest * 0 + agbstd_3).reshape(len(sampling_axis_north), len(sampling_axis_east)),
            ]

            upper_left_easting_coord = sampling_axis_east[0]  # i.e. horizontal
            sampling_step_east_west = sampling_axis_east[1] - sampling_axis_east[0]
            upper_left_northing_coord = sampling_axis_north[0]  # i.e. vertical
            sampling_step_north_south = sampling_axis_north[1] - sampling_axis_north[0]
            geotransform = [
                upper_left_easting_coord,
                sampling_step_east_west,
                0,
                upper_left_northing_coord,
                0,
                sampling_step_north_south,
            ]

            lut_cal_curr_row = [sampling_axis_east[0], sampling_axis_east[-1], sampling_axis_north[-1], sampling_axis_north[0], 1]
            self.lut_cal = np.vstack([self.lut_cal, lut_cal_curr_row])

            cal_out_intermediate = os.path.join(
                temp_agb_folder,
                'cal_estimated_from_block_{}'.format(current_block_index),
            )

            curr_cal_path = tiff_formatter(
                roi_vals,
                cal_out_intermediate,
                geotransform,
                multi_layers_tiff=True,
                gdal_data_format=gdal.GDT_Float32,
            )

            self.lut_cal_paths.append(curr_cal_path)

            sigma_tab = np.load(
                os.path.join(
                    sigma_tab_pixel_folder,
                    'sigma_tab_pixels_block_{}.npy'.format(current_block_index),
                )
            )
            theta_tab = np.load(
                os.path.join(
                    theta_tab_pixel_folder,
                    'theta_tab_pixels_block_{}.npy'.format(current_block_index),
                )
            )
            acqid_tab = np.load(
                os.path.join(
                    acq_tab_pixel_folder, 'acq_tab_pixels_block_{}.npy'.format(current_block_index)
                )
            )
            polid_tab = np.load(
                os.path.join(
                    polid_tab_pixel_folder,
                    'polid_tab_pixels_block_{}.npy'.format(current_block_index),
                )
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
                        (
                            forward(sigma_tab)
                            - par_tabs[0]
                            - par_tabs[2] * forward(np.cos(theta_tab))
                        )
                        * par_tabs[1],
                        axis=0,
                    )
                    / np.nansum(par_tabs[1] ** 2, axis=0)
                )
            )

            # full resolution values (estimate and stds)
            agbstd_2 = agbstd_2.reshape(len(sampling_axis_north), len(sampling_axis_east))
            agbstd_3 = agbstd_3.reshape(len(sampling_axis_north), len(sampling_axis_east))

            pix_vals = [
                agbest_tab.reshape(len(pixel_axis_north), len(pixel_axis_east)),
                (agbest_tab * 0 + agbstd_1).reshape(len(pixel_axis_north), len(pixel_axis_east)),
                sp.interpolate.griddata(
                    (sampling_mesh_east[~np.isnan(agbstd_2)], sampling_mesh_north[~np.isnan(agbstd_2)]),
                    agbstd_2[~np.isnan(agbstd_2)],
                    (pixel_mesh_east, pixel_mesh_north),
                    method='nearest',
                )
                * (~np.isnan(agbest_tab.reshape(len(pixel_axis_north), len(pixel_axis_east)))),
                sp.interpolate.griddata(
                    (sampling_mesh_east[~np.isnan(agbstd_3)], sampling_mesh_north[~np.isnan(agbstd_3)]),
                    agbstd_3[~np.isnan(agbstd_3)],
                    (pixel_mesh_east, pixel_mesh_north),
                    method='nearest',
                )
                * (~np.isnan(agbest_tab.reshape(len(pixel_axis_north), len(pixel_axis_east)))),
            ]

            quality_layer = np.sqrt(pix_vals[1] ** 2 + pix_vals[2] ** 2 + pix_vals[3] ** 2)

            agb_est_ground_out_file_curr = os.path.join(
                temp_agb_folder, 'agb_est_block_{}.tif'.format(current_block_index)
            )
            geotransform_curr = [pixel_axis_east[0], pixel_size_east, 0, pixel_axis_north[0], 0, pixel_size_north]

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
                lon_min, lat_min = getattr(self.e7g_intermediate, subgrid_code).xy2lonlat(
                    min(pixel_axis_east), min(pixel_axis_north)
                )
                lon_max, lat_max = getattr(self.e7g_intermediate, subgrid_code).xy2lonlat(
                    max(pixel_axis_east), max(pixel_axis_north)
                )
            except Exception as e:
                logging.error(
                    'Cannot recognize input FNF Equi7 mask "{}" sub-grid folder name :'.format(
                        subgrid_code
                    )
                    + str(e),
                    exc_info=True,
                )
                raise

            bbox = [(lon_min, lat_min), (lon_max, lat_max)]
            ftiles = e7g_product.search_tiles_in_roi(bbox=bbox)

            equi7_agb_est_outdir_curr = os.path.join(
                temp_agb_folder, 'eq7_agb_est_block_{}'.format(current_block_index)
            )
            os.makedirs(equi7_agb_est_outdir_curr)
            equi7_agb_est_outdir_temp = image2equi7grid(
                e7g_product,
                agb_est_ground_out_file_curr,
                equi7_agb_est_outdir_curr,
                gdal_path=self.gdal_path,
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
            block_finished_flag[current_block_index] = True
            # remove current par block from list
            block_order = block_order[block_order != current_block_index]
            # select next par block
            if len(block_order) > 0:
                # if there is at least one left, just take the next closest to CALdata
                current_block_index = block_order[0]
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

        for (
            tile_name,
            equi7_agb_est_out_tif_names_curr_tile,
        ) in equi7_agb_est_out_tif_names_per_tile.items():

            driver = gdal.Open(equi7_agb_est_out_tif_names_curr_tile[0], GA_ReadOnly)
            data_merged = np.nan * np.zeros((driver.RasterYSize, driver.RasterXSize))
            quality_merged = np.nan * np.zeros((driver.RasterYSize, driver.RasterXSize))
            driver = None

            for idx_tile, equi7_agb_est_out_tif_name_curr_tile in enumerate(
                equi7_agb_est_out_tif_names_curr_tile
            ):

                # merging current tile blocks togeter:
                driver = gdal.Open(equi7_agb_est_out_tif_name_curr_tile, GA_ReadOnly)

                if idx_tile == 0:
                    geotransform_out = driver.GetGeoTransform()
                elif geotransform_out != driver.GetGeoTransform():
                    err_str = 'Same equi7 tiles cannot have different geotrasform'
                    logging.error(err_str)
                    raise ValueError(err_str)

                data_merged = np.nanmean(
                    np.dstack((data_merged, driver.GetRasterBand(1).ReadAsArray())), axis=2
                )
                quality_merged = np.sqrt(
                    np.nanmean(
                        np.dstack(
                            (quality_merged ** 2, driver.GetRasterBand(2).ReadAsArray() ** 2)
                        ),
                        axis=2,
                    )
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
                global_agb_folder,
                sub_grid_folder_name,
                tile_name,
                'agb_est_' + tile_name + sub_grid_folder_name + '.tif',
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
            logging.info(
                '    equi7 tile blocks merged and savad to file{}'.format(output_folder_curr_tile)
            )

        if proc_conf.delete_temporary_files:
            try:
                shutil.rmtree(temp_proc_folder)
            except:
                pass
        logging.info('AGB: estimation ended correctly.\n')
