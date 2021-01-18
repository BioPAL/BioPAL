# this path is temporary:
from pathlib import Path
path_aux = Path( 'C:/bio/aux_for_agb_dev' )
      
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
    parse_agb_configuration_file,
    parse_coreprocessing_agb_configuration_file,
    write_coreprocessing_agb_configuration_file,
)

from biopal.io.data_io import (
    tiff_formatter,
)
from biopal.screen_calibration.screen_calibration import apply_calibration_screens
from biopal.geometry.utility_geometry import compute_and_oversample_geometry_auxiliaries
from biopal.ground_cancellation.ground_cancellation import ground_cancellation


# right now, let's do star import 

from biopal.agb.estimating_AGB import (
    transform_function,
    match_string_lists,
    sample_and_tabulate_data,
    fit_formula_to_random_subsets,
    save_human_readable_table,
    read_and_organise_3d_data,
    subset_iterable,
    map_space_variant_parameters,
    check_block_for_data_and_cal,
    compute_block_processing_order,
)
# %%
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
        stack_based_processing_obj = StackBasedProcessingAGB(
            self.configuration_file_xml,
            self.geographic_boundaries,
            self.geographic_boundaries_per_stack,
            self.gdal_path,
        )

        # Run Main APP #1: Stack Based Processing
        (coreprocessing_configuration_file_xml, lut_cal, lut_fnf, lut_stacks, equi7_initialization) = stack_based_processing_obj.run(
            input_file_xml
        )

        # Main APP #2: AGB Core Processing
        agb_processing_obj = CoreProcessingAGB(
            coreprocessing_configuration_file_xml,
            self.geographic_boundaries,
            self.gdal_path,
            lut_cal,
            lut_fnf,
            lut_stacks,
            equi7_initialization,
        )

        # Run Main APP #2: AGB Core Processing
        agb_processing_obj.run(input_file_xml)


class StackBasedProcessingAGB(Task):
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
        proc_inputs = parse_chains_input_file(input_file_xml)
        proc_conf = parse_agb_configuration_file(self.configuration_file_xml)

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
        coreprocessing_configuration_file_xml = r"C:\ARESYS_PROJ\BioPAL\biopal\conf\ConfigurationFile_CoreProcessingAGB_Mio.xml"
        
        return (
            coreprocessing_configuration_file_xml,
            LookupTableAGB(paths=lut_cal_paths, boundaries=lut_cal, progressive=None),
            LookupTableAGB(paths=lut_fnf_paths, boundaries=lut_fnf, progressive=None),
            LookupTableAGB(
                paths=lut_stacks_paths,
                boundaries=lut_stacks_boundaries,
                progressive=lut_progressive_stacks,
            ),
            equi7_initialization,
        )


class CoreProcessingAGB(Task):
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

        # AGB: Reading chains configuration files
        logging.info('AGB: Reading chains configuration files')
        check_if_path_exists(self.configuration_file_xml, 'FILE')
        proc_inputs = parse_chains_input_file(input_file_xml)
        proc_conf = parse_coreprocessing_agb_configuration_file(self.configuration_file_xml, proc_inputs.output_folder)

        # setting up directories and making sure that preprocessing has been run
        products_folder = os.path.join(proc_inputs.output_folder, 'Products')
        temp_proc_folder = os.path.join(products_folder, 'temporary_processing')
        if not (os.path.exists(temp_proc_folder)):
            error_message = '"temporary_processing" folder is not present in output: StackBasedProcessingAGB APP should be launched before CoreProcessingAGB '
            logging.error(error_message)
            raise RuntimeError(error_message)
        # check auxiliaries (equi7 initialization) and if not present, compute them
        self.check_auxiliaries(proc_conf)

        # preparing folders for final and temporary output
        global_agb_folder = os.path.join(products_folder, 'global_AGB')
        temp_agb_folder = os.path.join(temp_proc_folder, 'agb_estimation')
        os.makedirs(global_agb_folder)
        os.makedirs(temp_agb_folder)
        
        
        
        
        
        ### THIS PART READS THE SETUP
        
        '''
        
        Note: the xml file with setup and the associated reading function must
        be updated so that the formula, observables, and parameters all are defined
        
        
        '''
        proc_conf = parse_coreprocessing_agb_configuration_file(self.configuration_file_xml)
        # read and initialize all the parameters needed for the inversion
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
            sample_size_east,
            sample_size_north,
            block_spacing_east,
            block_spacing_north,
            block_size_east,
            block_size_north,
            number_of_subsets,
            geographic_grid_sampling,
            sub_grid_string,
        ) = initialize_inversion_parameters(
            self.equi7_sampling_intermediate,
            proc_inputs.geographic_grid_sampling,
            self.geographic_boundaries,
            proc_conf.AGB,
        )     

        #### Residual function defined as a list of function strings
        # (move this out of this code and into the xml file)
        
        # formula = ['l_hh + a_hh * agb_1_db + n_hh * cos_local_db + neg_sigma0_hh_db',
        #                         'l_hv + a_hv * agb_1_db + n_hv * cos_local_db + neg_sigma0_hv_db',
        #                         'l_vv + a_vv * agb_1_db + n_vv * cos_local_db + neg_sigma0_vv_db']
        
        # # experimental: testing the use of tomosar products
        # formula = ['l_hh + a_hh * agb_1_db + cos_tomo_theta_db + 10*np.log10(1-np.exp(-k_hh*tomo_h/cos_tomo2_theta)) + neg_cov_hh_30_db',
        #         'l_hv + a_hv * agb_1_db + cos_tomo_theta_db + 10*np.log10(1-np.exp(-k_hv*tomo_h/cos_tomo2_theta)) + neg_cov_hv_30_db',
        #         'l_vv + a_vv * agb_1_db + cos_tomo_theta_db + 10*np.log10(1-np.exp(-k_vv*tomo_h/cos_tomo2_theta)) + neg_cov_vv_30_db']
        
        
        #### Residual function defined as a list of function strings
        formula = proc_conf.AGB.residual_function.formula
        
        #### parameter names matching the unknown elements in the formula
        parameter_names = proc_conf.AGB.residual_function.formula_parameters.name
        
        #### parameter limits
        parameter_limits = proc_conf.AGB.residual_function.formula_parameters.limits
        
        #### flag indicating whether the parameters should be saved as a map
        parameter_save_as_map = proc_conf.AGB.residual_function.formula_parameters.save_as_map
        
        #### parameter variabilities 
        #  Defined over eight different dimensions (in this particular order):
        #   spatial - between samples
        #   forest class - between forest classes
        #   stack - between stacks, independent from where they come (equivalent to flagging all remaining flags)
        #   global cycle - between stacks if they come from different global cycle (equivalent to temporal of multiple months)
        #   heading - between stacks if they come from different headings
        #   swath - between stacks if they come from different 150-km swaths (equivalent to temporal of several weeks)
        #   subswath - between stacks if they come from different 50-km subswaths (equivalent to temporal of several days)
        #   azimuth - betwen stacks if they come from different 150-km images overlapping in azimuth (equivalent to temporal of a few seconds)
        number_of_parameters = len(proc_conf.AGB.residual_function.formula_parameters.name)

        parameter_variabilities = []
        for parameter_index in np.arange(number_of_parameters):
            parameter_variabilities.append(
            np.array([proc_conf.AGB.residual_function.formula_parameters.ChangesAcrossSamples[parameter_index],
                     proc_conf.AGB.residual_function.formula_parameters.ChangesAcrossForestClass[parameter_index],
                     proc_conf.AGB.residual_function.formula_parameters.ChangesAcrossStack[parameter_index],
                     proc_conf.AGB.residual_function.formula_parameters.ChangesAcrossGlobalCycle[parameter_index],
                     proc_conf.AGB.residual_function.formula_parameters.ChangesAcrossHeading[parameter_index],
                     proc_conf.AGB.residual_function.formula_parameters.ChangesAcrossSwath[parameter_index],
                     proc_conf.AGB.residual_function.formula_parameters.ChangesAcrossSubswath[parameter_index],
                     proc_conf.AGB.residual_function.formula_parameters.ChangesAcrossAzimuth[parameter_index]
                     ]
                     )
            )

                        
        #### parameter names matching the unknown elements in the formula
        #   these are the elements of the formula that will be estimated
        # parameter_names = ['l_hh','l_hv','l_vv',
        #                    'a_hh','a_hv','a_vv',
        #                    'n_hh','n_hv','n_vv',
        #                    'agb_1_db','agb_2_db','agb_3_db','agb_4_db',
        #                    'k_hh','k_hv','k_vv']
        
        
        #### parameter limits
        # for now, simply re-organise existing input
        # parameter_limits = [parameter_limits_l,parameter_limits_l,parameter_limits_l,
        #                     parameter_limits_a,parameter_limits_a,parameter_limits_a,
        #                     parameter_limits_n,parameter_limits_n,parameter_limits_n,
        #                     parameter_limits_w,[0,10],[0,10],[0,10],
        #                     [0,30*0.033],[0,30*0.033],[0,30*0.033]]
        
        #### flag indicating whether the parameters should be saved as a map
        # parameter_save_as_map = [False,False,False,
        #                          False,False,False,
        #                          False,False,False,
        #                          True,True,True,True,
        #                          False,False,False]
        
        #### parameter variabilities 
        #  Defined over eight different dimensions (in this particular order):
        #   spatial - between samples
        #   forest class - between forest classes
        #   stack - between stacks, independent from where they come (equivalent to flagging all remaining flags)
        #   global cycle - between stacks if they come from different global cycle (equivalent to temporal of multiple months)
        #   heading - between stacks if they come from different headings
        #   swath - between stacks if they come from different 150-km swaths (equivalent to temporal of several weeks)
        #   subswath - between stacks if they come from different 50-km subswaths (equivalent to temporal of several days)
        #   azimuth - betwen stacks if they come from different 150-km images overlapping in azimuth (equivalent to temporal of a few seconds)
        # for now, simply re-organise existing input
        # parameter_variabilities_orig = np.column_stack((np.ones((3,2))==0,np.array(parameter_variabilities_orig)[:,1:]))
        # parameter_variabilities_orig = parameter_variabilities_orig[np.array([0,0,0,1,1,1,2,2,2]),:]
        # parameter_variabilities = []
        # for parameter_variability in parameter_variabilities_orig:
        #     parameter_variabilities.append(parameter_variability)
        # parameter_variabilities.append(np.array([True,False,False,False,False,False,False,False]))
        # parameter_variabilities.append(np.array([True,False,False,False,False,False,False,False]))
        # parameter_variabilities.append(np.array([True,False,False,False,False,False,False,False]))
        # parameter_variabilities.append(np.array([True,False,False,False,False,False,False,False]))
        # parameter_variabilities.append(np.array([False,False,False,False,False,False,False,False]))
        # parameter_variabilities.append(np.array([False,False,False,False,False,False,False,False]))
        # parameter_variabilities.append(np.array([False,False,False,False,False,False,False,False]))
        
        #### Names for all elements of the formula that will be read from disk
        observable_names = proc_conf.AGB.residual_function.formula_observables.name
        
        #### flag indicating whether these observables are compulsory
        # if so, samples without this observable are rejected
        # (order of elements must match that in observable_names)
        observable_is_required = proc_conf.AGB.residual_function.formula_observables.is_required
        
        #### list with sources for observable data
        observable_sources = proc_conf.AGB.residual_function.formula_observables.source
        
        #### permissible intervals for source data (before transformation)
        observable_ranges = proc_conf.AGB.residual_function.formula_observables.ranges
        
        #### data transforms (here: tags indicating one of predefined functions)
        observable_transforms = proc_conf.AGB.residual_function.formula_observables.transform
        
        ### averaging method within current sampling area and stack
        #   although oftentimes 'mean' is used, this may not always be the case 
        #   e.g., for slope aspect angle v, arctan(mean(sin(v))/mean(cos(v))) 
        #   is a better averaging method as it is not as susceptible to 2pi ambiguities
        observable_averaging_methods = proc_conf.AGB.residual_function.formula_observables.averaging_method
       
        #### Names for all elements of the formula that will be read from disk
        # observable_names = ['neg_sigma0_hh_db','neg_sigma0_hv_db','neg_sigma0_vv_db',
        #                     'cos_local_db',
        #                     'agb_1_db','agb_2_db','agb_3_db','agb_4_db',
        #                     'lidar_h_db','tomo_h',
        #                     'neg_cov_hh_30_db','neg_cov_hv_30_db','neg_cov_vv_30_db',
        #                     'neg_cov_hh_40_db','neg_cov_hv_40_db','neg_cov_vv_40_db',
        #                     'cos_tomo_theta_db','ref_lidar_agb_db','cos_tomo2_theta']
        # note that agb is special: it is defined as both "observable" and "parameter", so it will be 
        # both read from external data and estimated 
        
        
        
        #### flag indicating whether these observables are compulsory
        # if so, samples without this observable are rejected
        # (order of elements must match that in observable_names)
        # observable_is_required = [True,True,True,
        #                           True,
        #                           False,False,False,False,
        #                           True,True,
        #                           True,True,True,
        #                           True,True,True,
        #                           True,True,True]
        
        #### list with sources for observable data
        # this part connects the observables to actual data
        # for now, it is assumed to consist of nested lists of two elements: 
        # path to a tiff file and band number
        # this may need to beorganised to make the code more user-friendly
        
        # current structure:
        #  observable_sources = [observable_1_sources,...,observable_n_sources,...]
        #   observable_n_sources = [stack_1_sources,...,stack_m_sources,...]
        #    stack_m_sources = [source_1,...,source_x,...]
        #     source_x = [tiff_file_path,band_id]
        
        # for now, simply re-organise existing input
        # observable_sources = []
        # for current_observable in ['sigma0_hh','sigma0_vh','sigma0_vv','theta']:
        #     stack_sources = []
        #     for current_stack in range(len(self.lut_stacks_paths)):
        #         file_sources = []
        #         current_path =  os.path.join(self.lut_stacks_paths[current_stack], current_observable) 
        #         for current_grid in os.listdir(current_path):
        #             current_path1 = os.path.join(current_path,current_grid)
        #             for current_tile in os.listdir(current_path1):
        #                 current_path2 = os.path.join(current_path1,current_tile)
        #                 for current_file in os.listdir(current_path2):
        #                     file_sources.append([os.path.join(current_path2,current_file),0])
        #         stack_sources.append(file_sources)
        #     observable_sources.append(stack_sources)
    
        # for current_observable_idx in range(4):
        #     file_sources = []
        #     for current_file in self.lut_cal_paths:
        #         file_sources.append([current_file,current_observable_idx])
        #     observable_sources.append([file_sources])
        # observable_sources = observable_sources+[
        #     [[[ str( path_aux.joinpath('lope_lidar','lidar_chm','EQUI7_AF050M','E045N048T3','lidar_chm_AF050M_E045N048T3.tif')),0]]],
        #     [[[ str( path_aux.joinpath('out_tomo_fh','BIOMASS_L2_20201118T123101','TOMO_FH','Products','global_FH','EQUI7_AF050M','E045N048T3','FH_EQUI7_AF050M_E045N048T3.tif')),0]]],
        #     [[[ str( path_aux.joinpath('out_tomo_GGS_50m_RES_200m','slice_30_m','EQUI7_AF050M','E045N048T3','cov_hh_30_m_EQUI7_AF050M_E045N048T3.tif')),0]]],
        #     [[[ str( path_aux.joinpath('out_tomo_GGS_50m_RES_200m','slice_30_m','EQUI7_AF050M','E045N048T3','cov_vh_30_m_EQUI7_AF050M_E045N048T3.tif')),0]]],
        #     [[[ str( path_aux.joinpath('out_tomo_GGS_50m_RES_200m','slice_30_m','EQUI7_AF050M','E045N048T3','cov_vh_30_m_EQUI7_AF050M_E045N048T3.tif')),0]]],
        #     [[[ str( path_aux.joinpath('out_tomo_GGS_50m_RES_200m','slice_40_m','EQUI7_AF050M','E045N048T3','cov_hh_40_m_EQUI7_AF050M_E045N048T3.tif')),0]]],
        #     [[[ str( path_aux.joinpath('out_tomo_GGS_50m_RES_200m','slice_40_m','EQUI7_AF050M','E045N048T3','cov_vh_40_m_EQUI7_AF050M_E045N048T3.tif')),0]]],
        #     [[[ str( path_aux.joinpath('out_tomo_GGS_50m_RES_200m','slice_40_m','EQUI7_AF050M','E045N048T3','cov_vh_40_m_EQUI7_AF050M_E045N048T3.tif')),0]]],
        #     [[[ str( path_aux.joinpath('out_tomo_GGS_50m_RES_200m','theta','EQUI7_AF050M','E045N048T3','theta_AF050M_E045N048T3.tif')),0]]],
        #     [[[ str( path_aux.joinpath('lope_lidar','lidar_agb','EQUI7_AF050M','E045N048T3','lidar_AGB_AF050M_E045N048T3.tif')),0]]],
        #     [[[ str( path_aux.joinpath('out_tomo_GGS_50m_RES_200m','theta','EQUI7_AF050M','E045N048T3','theta_AF050M_E045N048T3.tif')),0]]],
        #     ]
        
        
        
        # #### permissible intervals for source data (before transformation)
        # observable_ranges = [[1e-8,20],[1e-8,20],[1e-8,20],
        #                      [20*np.pi/180,60*np.pi/180],
        #                      [1,700],[0,10],[0,10],[0,10],
        #                      [0.1,100],[0.1,100],
        #                      [1e-8,20],[1e-8,20],[1e-8,20],
        #                      [1e-8,20],[1e-8,20],[1e-8,20],
        #                      [20*np.pi/180,60*np.pi/180],
        #                      [1,700],
        #                      [20*np.pi/180,60*np.pi/180]]
        
        # #### data transforms (here: tags indicating one of predefined functions)
        # observable_transforms = ['-db','-2db','-db',
        #                          'cosdb',
        #                          'db','none','none','none',
        #                          'db','none',
        #                          '-db','-2db','-db',
        #                          '-db','-2db','-db',
        #                          'cosdb',
        #                          'db','cos']
        
        # ### averaging method within current sampling area and stack
        # #   although oftentimes 'mean' is used, this may not always be the case 
        # #   e.g., for slope aspect angle v, arctan(mean(sin(v))/mean(cos(v))) 
        # #   is a better averaging method as it is not as susceptible to 2pi ambiguities
        # observable_averaging_methods = ['mean','mean','mean',
        #                                 'mean',
        #                                 'mean','mean','mean','mean',
        #                                 'mean','mean',
        #                                 'mean','mean','mean',
        #                                 'mean','mean','mean',
        #                                 'mean','mean','mean'] 
        
        
        
        
        

        ### read additional polygons (e.q., cals not on a grid)
        # (read from xml file or external shapefiles?)
        additional_sampling_polygons = []
        
        
        
        ### DERIVED QUANTITIES
        
        # some vectors flagging different types of observables
        observables_in_formula = np.any(match_string_lists(formula,observable_names)>=0,axis=0)
        observables_in_parameters = np.any(match_string_lists(parameter_names,observable_names)>=0,axis=0)
        observables_in_formula_or_parameters = observables_in_formula | observables_in_parameters
        
        
        
        # link parameters to observables in case the same name is observed in both
        positions_in_observable_list,positions_in_parameter_list = np.where(match_string_lists(observable_names,parameter_names)>=0)
        parameter_position_in_observable_vector = -1 * np.ones(len(parameter_names),dtype='int32')
        parameter_position_in_observable_vector[positions_in_parameter_list] = positions_in_observable_list
        
        # lists with paths to files that will be merged
        parameter_map_pathlists = [[],[],[],
                                   [],[],[],
                                   [],[],[],
                                   [],[],[],[],
                                   [],[],[]]
        
        # ### SOME MANIPULATIONS TO SPEED  UP THE INVERSION
        # # (remove in final version)
        
        # # take out the observables and parameters that are within formula or within the defined observables
        # # (no need for now to read unnecessary data)
        # observable_names = subset_iterable(observable_names,observables_in_formula_or_parameters,False)
        # observable_averaging_methods = subset_iterable(observable_averaging_methods,observables_in_formula_or_parameters,False)
        # observable_transforms = subset_iterable(observable_transforms,observables_in_formula_or_parameters,False)
        # observable_ranges = subset_iterable(observable_ranges,observables_in_formula_or_parameters,False)
        # observable_sources = subset_iterable(observable_sources,observables_in_formula_or_parameters,False)
        # observable_is_required = subset_iterable(observable_is_required,observables_in_formula_or_parameters,False)
        
        # # change the number of subsets
        # # (with the current setting, the subset size is 100% anyway, so no point running 5 subsets)
        # number_of_subsets = 1
        
        # remove in final versio
        fraction_of_samples = 0.8
        
        
        
        
        
        
        # checking that none of observable and parameter names can be contained within each other
        # (parsing of the formula will fail if this is the case)
        all_names = parameter_names + observable_names
        for current_name in all_names:
            for current_name2 in all_names:
                if ((current_name!=current_name2) & \
                    ((current_name.find(current_name2)>=0) | \
                     (current_name2.find(current_name)>=0))):
                    logging.error('AGB: parameter and observable names must be unique and cannot contain each other ({} and {} do not fulfill this requirement).'.format(current_name,current_name2))
            
                    
        
       
            

        ### PREPARING STACK INFO
        # read acquisition info table
        stack_info_table = self.lut_progressive_stacks
        stack_info_table_columns = ['stack_id','global_cycle_id','heading_degrees','swath_id','subswath_id','azimuth_id']
        
       
        ### PREPARING PARAMETER GRID
        # parameter block mesh and flattened coordinate vectors
        block_mesh_east, block_mesh_north = np.meshgrid(np.arange(first_pixel_east, last_pixel_east, block_spacing_east), np.arange(first_pixel_north, last_pixel_north, block_spacing_north))
        block_corner_coordinates_east, block_corner_coordinates_north = [
            np.float64(x.flatten()) for x in [block_mesh_east, block_mesh_north]
        ]
        # count the blocks and set up a "finished" flag vector
        number_of_blocks = len(block_corner_coordinates_east.flatten())
        
        
        # Compute the starting block and the block processing order
        block_order = compute_block_processing_order(
                block_corner_coordinates_east,
                block_corner_coordinates_north, 
                block_size_east, 
                block_size_north,
                self.lut_cal, # right now, this uses boundaries but in the future it should be capable of using polygons (including additional_sampling_polygons)
                self.lut_stacks_boundaries # right now, this uses boundaries but in the future it should be capable of using polygons
            )

        
        
        
        
        ### RUNNING PARAMETER BLOCKS
        # block status vector
        block_status = np.zeros(number_of_blocks)
        # flag explanation
        # -1: error (see log for more details)
        # 0: block not run
        # 1: block skipped due to insufficient data
        # 2: block skipped due to lack of sufficient subsets
        # 3: block successfully finished
        
        
        
        for counter_blocks_run,current_block_index in enumerate(block_order):
        
            # flag indicating whether to skip the current block or not
            # (must be reset upon each iteration)
            skip_current_block = False
            
            logging.info('AGB: Running block {} out of {} (block ID: {})'.format(counter_blocks_run+1, number_of_blocks, current_block_index))
            
       
            # %% ### CREATING SAMPLING GRID AND TESTING FOR DATA
            try:
                logging.info('AGB: Creating sampling grid and checking for data.')
            
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
                
                
                # calculate the total number of samples
                number_of_samples = number_of_samples_on_grid + len(additional_sampling_polygons)
                
                # checking the number of samples
                if number_of_samples < proc_conf.AGB.min_number_of_rois:
                    logging.info(
                        '... skipping block #{} because the number of samples #{} cannot be less than #{}'.format(
                            current_block_index, number_of_samples_on_grid, proc_conf.AGB.min_number_of_rois))
                    skip_current_block = True
                    block_status[counter_blocks_run] = 1
                    
                else:
                    # examine stack and calibration data
                    (
                        block_has_data,
                        block_has_cal
                        ) = check_block_for_data_and_cal(
                            current_block_extents,
                            self.lut_stacks_boundaries,
                            self.lut_cal)
                            
                    
                    # checking availability of calibration and stack data        
                    if ~np.any(block_has_data) or ~np.any(block_has_cal):
                        logging.info('... skipping block #{} due to lack of valid data or calibration points'.format(
                                current_block_index))
                        skip_current_block = True
                        block_status[counter_blocks_run] = 1
                        
                    
            except Exception as e:
                logging.error('AGB: error during sampling grid preparation.' + str(e), exc_info=True)
                block_status[counter_blocks_run] = -1
                raise
                
            if skip_current_block:
                continue
            # %% ### TABULATING DATA
            try:
                logging.info('AGB: tabulating data...')
                
                # pixel axes for current block
                pixel_axis_east = np.arange(current_block_extents[0], current_block_extents[1], pixel_size_east)
                pixel_axis_north = np.arange(current_block_extents[2], current_block_extents[3], pixel_size_north)
                
                
                # output equi7 projection info
                #   here, we assume that the output tile and subtile will be that of the first observable source
                #   that is covered by the current block
                equi7_info_source_path = observable_sources[0][np.where(block_has_data)[0][0]][0][0]
                equi7_subtile_name,equi7_tile_name = [x.split('.')[0] for x in equi7_info_source_path.split('_')[-2:]]
                equi7_subgrid_code = equi7_subtile_name[:2]
                equi7_projection_string = get_projection_from_path(equi7_info_source_path)
                current_geotransform = [pixel_axis_east[0], pixel_size_east, 0, pixel_axis_north[0], 0, pixel_size_north]        
                equi7_product = Equi7Grid(geographic_grid_sampling)
                
                
                
                
                # tabulation
                (
                    observable_table,
                    observable_names,
                    identifier_table,
                    identifier_names,
                    parameter_position_table,
                    parameter_position_table_columns,
                    parameter_tables,
                    parameter_table_columns,
                ) = sample_and_tabulate_data(
                        current_block_extents, # extent of the current area for which the table is created
                        pixel_axis_east, # east north axes onto which data are interpolated
                        pixel_axis_north,
                        sampling_axis_east,
                        sampling_axis_north,
                        sample_size_east, # east north extents of the samples
                        sample_size_north,
                        additional_sampling_polygons, # additional arbitrarily shaped polygons
                        block_has_data, # flags whether each stack is in current block
                        stack_info_table, # info table with stack properties (stack id, headings, etc., defining the acquisition parameters)
                        stack_info_table_columns, # column names for the abovementioned table
                        self.lut_fnf,
                        [[x,0] for x in self.lut_fnf_paths],
                        observable_names, # observable names in formula
                        observable_is_required,
                        observable_sources, # paths and band ids
                        observable_transforms, # transform function to apply to observable
                        observable_averaging_methods, # averaging method (most commonly 'mean', but could be other if required (e.g., for slope aspect angle))
                        observable_ranges, # permitted ranges, outside those the observable is set to nan
                        parameter_names, # parameter names in formula
                        parameter_limits, # permissible parameter intervals
                        parameter_variabilities, # parameter variabilities across all dimensions
                        number_of_subsets, # number of subsets to use (used to allocate columns in parameter tables)
                )
                
                # checking if the tables contain any data
                if observable_table.shape[0]==0:
                    logging.info('... skipping block #{} due to lack of valid data points'.format(
                            current_block_index))
                    skip_current_block = True
                    block_status[counter_blocks_run] = 1
                 
            except Exception as e:
                logging.error('AGB: error during data sampling and tabulation.' + str(e), exc_info=True)
                block_status[counter_blocks_run] = -1
                raise
                
            if skip_current_block:
                continue
                
            # %% ### FITTING MODEL TO TABULATED DATA  AND SAVNG
            
            try:
                    
                logging.info('AGB: fitting model to data...')
                
                (
                    parameter_tables,
                    space_invariant_parameter_table,
                    space_invariant_parameter_names,
                    space_variant_parameter_table,
                    space_variant_parameter_names,
                ) = fit_formula_to_random_subsets(
                        formula,
                        number_of_subsets,
                        observable_table,
                        observable_names,
                        identifier_table,
                        identifier_names,
                        parameter_position_table, 
                        parameter_names,
                        parameter_tables,
                        parameter_table_columns,
                        parameter_variabilities,
                        proc_conf.AGB.fraction_of_cal_per_test/100*fraction_of_samples,
                        proc_conf.AGB.fraction_of_roi_per_test/100*fraction_of_samples,
                        proc_conf.AGB.min_number_of_cals_per_test,
                        proc_conf.AGB.min_number_of_rois_per_test,
                )
                
                # checking if the tables contain any data
                if len(space_invariant_parameter_table)==0:
                    logging.info('... skipping block #{} due to no subsets found.'.format(
                            current_block_index))
                    skip_current_block = True
                    block_status[counter_blocks_run] = 2
                else:
                    logging.info('AGB: saving data to tables.')
                    # save data
                    parameter_property_names = parameter_table_columns[0][-number_of_subsets-4:]
                    line_number_string = ['row']
                    
                    # select formatting for the output tables
                    table_delimiter = '\t'
                    table_precision = 3
                    table_column_width = 25
                    data_type_lut = [[line_number_string,identifier_names,parameter_position_table_columns,parameter_property_names,observable_names,parameter_names],
                                     ['d','d','d','f','f','f']]
                    
                    # save observable table (includes indices to parameter tables, replicated if necessary)
                    curr_path = os.path.join(
                            temp_agb_folder,
                            'observable_table_block_{}.txt'.format(current_block_index))
                    curr_table = np.column_stack((
                        np.arange(observable_table.shape[0]),
                        identifier_table,
                        observable_table,
                        parameter_position_table))
                    curr_column_names = line_number_string+identifier_names+observable_names+parameter_position_table_columns
                    save_human_readable_table(curr_path,curr_table,curr_column_names,data_type_lut,table_delimiter,table_precision,table_column_width)                
                    
                    # save results table (with parameter estimates replicated across samples etc, if necessary; no indices connecting these to parameter tables)
                    curr_path = os.path.join(
                            temp_agb_folder,
                            'results_table_block_{}.txt'.format(current_block_index))
                    curr_table = np.column_stack((
                        np.arange(observable_table.shape[0]),
                        identifier_table,
                        observable_table,
                        space_invariant_parameter_table,
                        space_variant_parameter_table))
                    curr_column_names = line_number_string+identifier_names+observable_names+space_invariant_parameter_names+space_variant_parameter_names
                    save_human_readable_table(curr_path,curr_table,curr_column_names,data_type_lut,table_delimiter,table_precision,table_column_width)                
                    
                    # save parameter tables (minimalistic tables with only the most necessary info)
                    for parameter_idx,parameter_name in enumerate(parameter_names):
                        curr_path = os.path.join(
                            temp_agb_folder,
                            'parameter_{}_table_block_{}.txt'.format(parameter_name,current_block_index))
                        curr_table = np.column_stack((
                            np.arange(parameter_tables[parameter_idx].shape[0]),
                            parameter_tables[parameter_idx]))
                        curr_column_names = np.concatenate((np.array(line_number_string),parameter_table_columns[parameter_idx]))
                        save_human_readable_table(curr_path,curr_table,curr_column_names,data_type_lut,table_delimiter,table_precision,table_column_width)                

            except Exception as e:
                logging.error('AGB: error during parameter estimation or data saving.' + str(e), exc_info=True)
                block_status[counter_blocks_run] = -1
                raise
                
            
            if skip_current_block:
                continue
                
                
            # %% ### READING IMAGE DATA
            try:
                
                logging.info('AGB: reading data images (only necessary files).')
                
                # take out the observables that are in formula and not among space variant parameters
                observables_for_mapping = np.any(match_string_lists(formula,observable_names)>=0,axis=0) & \
                            ~np.any(match_string_lists(space_variant_parameter_names,observable_names)>=0,axis=0)
                
                (
                    forest_class_3d,
                    observables_3d,
                    observables_3d_names,
                    space_invariant_parameters_3d,
                    space_invariant_parameters_3d_names,
                    identifiers_3d,
                    identifiers_3d_names,
                ) = read_and_organise_3d_data(
                        current_block_extents,
                        block_has_data,
                        pixel_axis_north,
                        pixel_axis_east,
                        stack_info_table,
                        stack_info_table_columns,
                        subset_iterable(observable_names,observables_for_mapping,False),
                        subset_iterable(observable_sources,observables_for_mapping,False),
                        subset_iterable(observable_transforms,observables_for_mapping,False),
                        subset_iterable(observable_averaging_methods,observables_for_mapping,False),
                        subset_iterable(observable_ranges,observables_for_mapping,False),
                        subset_iterable(observable_is_required,observables_for_mapping,False),
                        [[x,0] for x in self.lut_fnf_paths],
                        self.lut_fnf,
                        identifier_table[:,2],
                        identifier_table[:,1],
                        space_invariant_parameter_table,
                        space_invariant_parameter_names,
                        mask_out_area_outside_block = True
                )
                    
                # skipping to next if all parameters and/or forest class data are missing
                if np.all(forest_class_3d==0) | np.all(
                        np.array(
                            [np.all(
                                np.isnan(
                                    curr_space_invariant_parameter_3d
                                    )
                                ) for curr_space_invariant_parameter_3d in space_invariant_parameters_3d]
                            )
                        ):
                    logging.info('... skipping block #{} due to no image data read.'.format(
                            current_block_index))
                    skip_current_block = True
                    block_status[counter_blocks_run] = 1
                    
                    
            except Exception as e:
                logging.error('AGB: error during image reading.' + str(e), exc_info=True)
                block_status[counter_blocks_run] = -1
                raise
                
                
            if skip_current_block:
                continue
                
            # %% ### ESTIMATING SPACE-VARIANT PARAMTERS
            
            try:
                          
                logging.info('AGB: creating space variant parameter images.')

                # take out the observables that are in formula and not among space variant parameters
                parameters_for_mapping = np.any(match_string_lists(formula,parameter_names)>=0,axis=0) & \
                            ~np.any(match_string_lists(space_invariant_parameter_names,parameter_names)>=0,axis=0)
                
                (
                    space_variant_parameters_3d,
                    space_variant_parameters_3d_names,
                ) = map_space_variant_parameters(
                        formula,
                        forest_class_3d,
                        observables_3d,
                        observables_3d_names,
                        space_invariant_parameters_3d,
                        space_invariant_parameters_3d_names,
                        identifiers_3d,
                        identifiers_3d_names,
                        subset_iterable(parameter_names,parameters_for_mapping,False),
                        subset_iterable(parameter_variabilities,parameters_for_mapping,False),
                        subset_iterable(parameter_limits,parameters_for_mapping,False))
                
                # skipping to next if all parameters and/or forest class data are missing
                if np.all(
                        np.array(
                            [np.all(
                                np.isnan(
                                    curr_space_variant_parameter_3d
                                    )
                                ) for curr_space_variant_parameter_3d in space_variant_parameters_3d]
                            )
                        ):
                    logging.info('... skipping block #{} due to no space variant parameter maps created'.format(
                            current_block_index))
                    skip_current_block = True
                    block_status[counter_blocks_run] = 1
                      
                   
            except Exception as e:
                logging.error('AGB: error during space variant parameter mapping.' + str(e), exc_info=True)
                block_status[counter_blocks_run] = -1
                raise
                
                
            if skip_current_block:
                continue
                
            # %% ### ESTIMATING OTHER PARAMETERS
            
            try:
                
                # this is where other things such as error parameters are estimated
                logging.info('AGB: estimating other parameters to be implemented... (includes estimation of error)')
                
                
                # for now, the three error elements are hard coded
                additional_parameters_3d = [space_variant_parameter_3d*0 for space_variant_parameter_3d in space_variant_parameters_3d*3]
                additional_parameters_3d_names = ['agb_2_db','agb_3_db','agb_4_db']
            
            except Exception as e:
                logging.error('AGB: error during estimation of other parameters.' + str(e), exc_info=True)
                block_status[counter_blocks_run] = -1
                raise
                
                
            if skip_current_block:
                continue
                
            # %% ### SAVING SPACE-VARIANT PARAMTERS:
                     
            try:         
                logging.info('AGB: saving selected parameters to images.')
                
                # combine all available mapped parameters
                all_parameters_3d_names = space_invariant_parameters_3d_names +\
                    space_variant_parameters_3d_names+\
                        additional_parameters_3d_names
                all_parameters_3d = space_invariant_parameters_3d +\
                    space_variant_parameters_3d +\
                        additional_parameters_3d
                        
                
                
                for parameter_idx,parameter_name in enumerate(parameter_names):
                    if parameter_save_as_map[parameter_idx]:
                        
                        current_image_to_write = all_parameters_3d[np.where(np.array(all_parameters_3d_names)==parameter_name)[0][0]]
                        
                        current_file_path = os.path.join(
                            temp_agb_folder,parameter_name)
                        
                        # if the parameter is also an observable, we need to apply inverse transform
                        current_position_in_observable_vector = parameter_position_in_observable_vector[parameter_idx]
                        if current_position_in_observable_vector!=-1:
                            current_image_to_write = transform_function(
                                current_image_to_write,
                                observable_ranges[current_position_in_observable_vector],
                                observable_transforms[current_position_in_observable_vector],
                                do_forward = False
                            )
                            current_file_path += '_backtransf'
                            
                        current_file_path +='_block_%d' % (current_block_index) + '.tif'# + '_'+equi7_subtile_name+'_'+equi7_tile_name
                        
                            
                            
                        # save tiff file
                        current_file_path = tiff_formatter(
                            [x for x in current_image_to_write.transpose([2,0,1])],
                            current_file_path,
                            current_geotransform,
                            gdal_data_format=gdal.GDT_Float32,
                            projection=equi7_projection_string,
                            multi_layers_tiff=True,
                        )
                        
                        
                        # get extents of the current EQUI7 tile             
                        # [(left, lower), (right, upper)]
                        try:
                            lon_min, lat_min = getattr(self.e7g_intermediate, equi7_subgrid_code).xy2lonlat(
                                min(pixel_axis_east), min(pixel_axis_north)
                            )
                            lon_max, lat_max = getattr(self.e7g_intermediate, equi7_subgrid_code).xy2lonlat(
                                max(pixel_axis_east), max(pixel_axis_north)
                            )
                        except Exception as e:
                            logging.error(
                                'Cannot recognize input FNF Equi7 mask "{}" sub-grid folder name :'.format(
                                    equi7_subgrid_code
                                )
                                + str(e),
                                exc_info=True,
                            )
                            raise
                        
                        output_equi7_file_path = image2equi7grid(
                            equi7_product,
                            current_file_path,
                            temp_agb_folder,
                            gdal_path=self.gdal_path,
                            ftiles=equi7_product.search_tiles_in_roi(bbox=[(lon_min, lat_min), (lon_max, lat_max)]),
                            accurate_boundary=False,
                            tile_nodata=np.nan,
                        )
                        
                    
                        if current_position_in_observable_vector!=-1:
                            observable_sources[current_position_in_observable_vector][0].append([current_file_path,0])
                        parameter_map_pathlists[parameter_idx].append([output_equi7_file_path[0],0])
                                
                                
                                
                        # self.lut_cal_paths.append(output_file_path)
                self.lut_cal = np.row_stack((self.lut_cal, np.concatenate((current_block_extents[np.array([0,1,3,2])],np.zeros(1)))))
                    
                        
                      
                
                skip_current_block = False
                block_status[counter_blocks_run] = 3
                
                
            except Exception as e:
                logging.error('AGB: error during saving of maps.' + str(e), exc_info=True)
                raise
                    
        
        # %% FINAL MERGING OF THE IMAGES
            
        
        
        
        try:
            logging.info('AGB: creating wall-to-wall maps...')
            for parameter_idx,parameter_name in enumerate(parameter_names):
                if parameter_save_as_map[parameter_idx]:
                    
                    tiles_to_save = {}
                    for current_source in parameter_map_pathlists[parameter_idx]:
                            subtile_name,tile_name = [x.split('.')[0] for x in current_source[0].split('_')[-2:]]
                        
                            if tile_name in tiles_to_save.keys():
                                tiles_to_save[tile_name].append(current_source)
                            else:
                                tiles_to_save[tile_name] = [current_source]
                    
                    for tile_name,tile_sources in tiles_to_save.items():
                    
                        driver = gdal.Open(tile_sources[0][0], GA_ReadOnly)
                        data_merged = np.nan * np.zeros((driver.RasterYSize, driver.RasterXSize))
                        driver = None
                    
                        for idx_tile, current_tile_source in enumerate(
                            tile_sources
                        ):
                    
                            # merging current tile blocks togeter:
                            driver = gdal.Open(current_tile_source[0], GA_ReadOnly)
                    
                            if idx_tile == 0:
                                geotransform_out = driver.GetGeoTransform()
                            elif geotransform_out != driver.GetGeoTransform():
                                err_str = 'Same equi7 tiles cannot have different geotrasform'
                                logging.error(err_str)
                                raise ValueError(err_str)
                    
                            data_merged = np.nanmean(
                                np.dstack((data_merged, driver.GetRasterBand(current_tile_source[1]+1).ReadAsArray())), axis=2
                            )
                            
                            driver = None
                    
                        
                        current_merged_file_path = os.path.join(
                            global_agb_folder,
                            equi7_subtile_name,
                            tile_name,
                            parameter_name,
                        )
                        
                        current_position_in_observable_vector = parameter_position_in_observable_vector[parameter_idx]
                        if current_position_in_observable_vector!=-1:
                            current_merged_file_path += '_backtransf_'
                        current_merged_file_path += tile_name + '_' + equi7_subtile_name + '.tif'
                        
                        
                        
                        if not os.path.exists(os.path.dirname(current_merged_file_path)):
                            os.makedirs(os.path.dirname(current_merged_file_path))
                            
                        current_merged_file_path = tiff_formatter(
                            [data_merged],
                            current_merged_file_path,
                            geotransform_out,
                            gdal_data_format=gdal.GDT_Float32,
                            projection=equi7_projection_string,
                            multi_layers_tiff=True,
                        )
                        logging.info(
                            "... successful for parameter '{}' and filename '{}'".format(parameter_name, current_merged_file_path)
                        )
                    
                    # if proc_conf.delete_temporary_files:
                    #     try:
                    #         shutil.rmtree(temp_proc_folder)
                    #     except:
                    #         pass
            logging.info('AGB: estimation ended correctly.\n')
        
        
        except Exception as e:
            logging.error('AGB: error during creation of wall-to-wall maps.' + str(e), exc_info=True)
            raise