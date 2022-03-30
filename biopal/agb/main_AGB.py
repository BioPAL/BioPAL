# SPDX-FileCopyrightText: BioPAL <biopal@esa.int>
# SPDX-License-Identifier: MIT

import os
import numpy as np
import logging
from shapely.geometry import Polygon
from osgeo.gdalconst import GA_ReadOnly
from osgeo import gdal
from equi7grid.equi7grid import Equi7Grid
from equi7grid.image2equi7grid import image2equi7grid
from scipy.signal import convolve2d

# biomassL2 processor imports
from biopal.agb.processing_AGB import (
    LookupTableAGB,
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
    start_logging,
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
    set_gdal_paths,
    radiometric_correction_beta_to_sigma,
)
from biopal.geocoding.geocoding import (
    geocoding,
    geocoding_init,
)
from biopal.io.xml_io import (
    parse_input_file,
    write_input_file,
    parse_configuration_file,
    write_configuration_file,
    core_processing_agb,
)

from biopal.io.data_io import tiff_formatter
from biopal.screen_calibration.screen_calibration import apply_calibration_screens
from biopal.geometry.utility_geometry import compute_and_oversample_geometry_auxiliaries
from biopal.ground_cancellation.ground_cancellation import ground_cancellation

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


class InvalidInputError(ValueError):
    pass


class AboveGroundBiomass(Task):
    """AGB main APP Above Ground Biomass
    
    run this APP to execute the complete Above Ground Biomass (AGB) processing chain.
    
    AboveGroundBiomass is composed by two sub APPS automatically called in sequence  
    when standard launch is performed:
    StackBasedProcessingAGB -> CoreProcessingAGB

    Refer to dataset_query, StackBasedProcessingAGB and CoreProcessingAGB documentation 
    and Tutorial section for step by step run.

    Attributes
    ----------
    configuration_file : str
        path of the Configuration_File.xml file

    Methods
    -------
    run( input_file_path )
        run the AboveGroundBiomass processing
    name : str
        name of the APP

    See Also
    --------
    biopal.dataset_query.dataset_query.dataset_query: it's the APP to be called before this APP
    StackBasedProcessingAGB : it's the first of the two sub-APPs called by AboveGroundBiomass preocessor
    CoreProcessingAGB : it's  the second of the two sub-APPs called by AboveGroundBiomass processor

    Examples
    --------
    Manual AGB chain execution

    >>> from biopal.dataset_query.dataset_query import dataset_query
    >>> from biopal.agb.main_AGB import AboveGroundBiomass
    >>> dq_obj = dataset_query()
    >>> input_file_up = dq_obj.run( input_file )
    >>> chain_obj = AboveGroundBiomass( configuration_file )
    >>> chain_obj.run( input_file_up )

    - input_file: path of the BioPAL input file
    - input_file_up: same of input_file with also the "stack_based_processing" section
    - configuration_file: path of the BioPAL configuration file
    """

    def __init__(
        self, configuration_file,
    ):
        super().__init__(configuration_file)

    def _run(self, input_file):

        # Main APP #1: Stack Based Processing
        stack_based_processing_obj = StackBasedProcessingAGB(self.configuration_file)

        # Run Main APP #1: Stack Based Processing
        input_file_updated, configuration_file_updated = stack_based_processing_obj.run(input_file)

        # Main APP #2: AGB Core Processing
        agb_processing_obj = CoreProcessingAGB(configuration_file_updated)

        # Run Main APP #2: AGB Core Processing
        agb_processing_obj.run(input_file_updated)


class StackBasedProcessingAGB(Task):
    """StackBasedProcessingAGB APP 
    
    StackBasedProcessingAGB APP is the first of the two sub-APPs called by AboveGroundBiomass processor.

    It performs stack-based operations to prepare inputs for the "CoreProcessingAGB" APP. 
    It is automatically called in the default run (see Getting Started  and Tutorial documentation section),
    or can be manually launched (stand alone) following this sequence:
    dataset_query -> StackBasedProcessingAGB -> CoreProcessingAGB

    Attributes
    ----------
    configuration_file : str
        path of the Configuration_File.xml file

    Methods
    -------
    run( input_file_path )
        run the StackBasedProcessingAGB APP
    name : str
        name of the APP

    See Also
    --------
    biopal.dataset_query.dataset_query.dataset_query : it's the first APP to be called in the manual sequence
    CoreProcessingAGB : it's the core APP that follows this APP in the call sequence

    Examples
    --------
    Manual AGB chain execution

    >>> from biopal.dataset_query.dataset_query import dataset_query
    >>> from biopal.AGB.main_AGB import StackBasedProcessingAGB, CoreProcessingAGB
    >>> dq_obj = dataset_query()
    >>> input_file_up1 = dq_obj.run( input_file )
    >>> sbp_obj = StackBasedProcessingAGB( config_file )
    >>> input_file_up2, config_file_up = sbp_obj.run( input_file_up1 )
    >>> agbp_obj = CoreProcessingAGB( config_file_up )
    >>> agbp_obj.run( input_file_up2 )

    - input_file: path of the BioPAL input file
    - input_file_up1: same of input_file with also the "stack_based_processing" section   
    - input_file_up2: same of input_file_up1 with also the "core_processing_agb" section
    - config_file: path of the BioPAL configuration file
    - config_file_up: same of config_file, with paths updated      
    """

    def __init__(self, configuration_file):
        super().__init__(configuration_file)

    def _run(self, input_file):

        # parse the input file and check for the needed sections in it:
        if isinstance(input_file, str):
            check_if_path_exists(input_file, "FILE")
            input_params_obj = parse_input_file(input_file)
        else:
            input_params_obj = input_file
        if (
            input_params_obj.L2_product is None
            or input_params_obj.output_specification is None
            or input_params_obj.dataset_query is None
            or input_params_obj.stack_based_processing is None
        ):
            error_message = [
                "StackBasedProcessingAGB APP input file should contain at least"
                ' "L2_product", "output_specification", "dataset_query" and "stack_based_processing" sections'
            ]
            raise InvalidInputError(error_message)

        # parse the configuration file and check for the needed sections in it:
        if isinstance(input_file, str):
            check_if_path_exists(self.configuration_file, "FILE")
            conf_params_obj = parse_configuration_file(self.configuration_file)
        else:
            conf_params_obj = self.configuration_file
        if (
            conf_params_obj.processing_flags is None
            or conf_params_obj.ground_cancellation is None
            or conf_params_obj.estimate_agb is None
        ):
            error_message = [
                "Configuration file for  StackBasedProcessingAGB APP should contain at least"
                ' "processing_flags", "ground_cancellation" and "estimate_agb" sections'
            ]
            raise InvalidInputError(error_message)

        # Initialize the logger, of the APP has been called in a stand-alone call
        if not hasattr(logging.getLoggerClass().root.handlers[0], "baseFilename"):
            start_logging(
                os.path.dirname(input_params_obj.output_specification.output_folder),
                input_params_obj.L2_product,
                "DEBUG",
                app_name="StackBasedProcessingAGB",
            )
        logging.info("AGB stack-based-processing APP starting\n")

        # get needed parameters from input and configuration files
        geographic_boundaries = input_params_obj.stack_based_processing.geographic_boundaries
        geographic_boundaries_per_stack = input_params_obj.stack_based_processing.geographic_boundaries_per_stack

        gdal_path, _ = set_gdal_paths(conf_params_obj.gdal.gdal_path, conf_params_obj.gdal.gdal_environment_path)

        ########################## INITIAL STEPS ##############################

        ### managing output folders:
        products_folder = os.path.join(input_params_obj.output_specification.output_folder, "Products")
        if conf_params_obj.processing_flags.save_breakpoints:
            breakpoints_output_folder = os.path.join(products_folder, "breakpoints")
            logging.info("AGB: Breakpoints will be saved into: " + breakpoints_output_folder)
            os.makedirs(breakpoints_output_folder)

        temp_output_folder = os.path.join(products_folder, "temp")
        logging.info("AGB: Temporary data folder:" + temp_output_folder + "\n")

        os.makedirs(temp_output_folder)

        ### get temporal date time of the input data (get the minimum date from all the stacks)
        time_tag_mjd_initial = get_min_time_stamp_repository(
            input_params_obj.dataset_query.L1C_repository, input_params_obj.stack_based_processing.stack_composition
        )

        # read the cals:
        cal_format = check_cal_format(input_params_obj.stack_based_processing.reference_agb_folder)
        if cal_format == "GeoJSON":
            cal_fnames = get_foss_cal_names(input_params_obj.stack_based_processing.reference_agb_folder)
            flag_cal = 0
            pass
        elif cal_format == "RASTER":
            flag_cal = 1
            cal_fnames = get_raster_cal_names(input_params_obj.stack_based_processing.reference_agb_folder)

            # LUT: CAL
            lut_cal_boundaries = np.zeros((len(cal_fnames), 5))
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

                lut_cal_boundaries[idx_cal, :] = [
                    east_min,
                    east_max,
                    min(north_in, north_out),
                    max(north_in, north_out),
                    flag_cal,
                ]
                lut_cal_paths.append(equi7_cal_fname)

                logging.info("using calibration file: " + equi7_cal_fname)

        ### initialize the equi7 sampling grid
        equi7_sampling_intermediate = choose_equi7_sampling(
            conf_params_obj.estimate_agb.intermediate_ground_averaging,
            conf_params_obj.estimate_agb.intermediate_ground_averaging / 2,
        )
        e7g_intermediate = Equi7Grid(equi7_sampling_intermediate)
        logging.info("EQUI7 Grid sampling used for intermediate products: {}".format(equi7_sampling_intermediate))

        # it is loaded if equi7, or converted to equi7 if TANDEM-X
        # it is a list containing all the loaded FNF-FTILES
        fnf_format = check_fnf_folder_format(input_params_obj.stack_based_processing.forest_mask_catalogue_folder)
        if fnf_format == "TANDEM-X":
            logging.info("Initial Forest mask is in TANDEM-X format, converting to equi7...")

            # conversion step 1: get the tiff names of current zone
            equi7_fnf_mask_fnames = fnf_tandemx_load_filter_equi7format(
                input_params_obj.stack_based_processing.forest_mask_catalogue_folder,
                e7g_intermediate,
                conf_params_obj.estimate_agb.intermediate_ground_averaging,
                temp_output_folder,
                gdal_path,
                geographic_boundaries,
                time_tag_mjd_initial,
            )

        elif fnf_format == "EQUI7":
            logging.info("Initial Forest mask reading and formatting...")

            # re format the input mask in equi7 with the output resolution:
            equi7_fnf_mask_fnames = fnf_equi7_load_filter_equi7format(
                input_params_obj.stack_based_processing.forest_mask_catalogue_folder,
                e7g_intermediate,
                conf_params_obj.estimate_agb.intermediate_ground_averaging,
                temp_output_folder,
                gdal_path,
            )

        # LUT: FNF
        lut_fnf_boundaries = np.zeros((len(equi7_fnf_mask_fnames), 4))
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

            lut_fnf_boundaries[fnf_idx, :] = [
                east_min,
                east_max,
                min(north_in, north_out),
                max(north_in, north_out),
            ]
            lut_fnf_paths.append(fnf_name_curr)

        # To be set from configuration file in the final version
        space_varying_ground_cancellation = False

        ########################## INITIAL STEPS END #############################

        ########################## STACK BASED STEPS ##############################

        # LookUp tables initialization
        number_of_stacks = len(input_params_obj.stack_based_processing.stack_composition.keys())
        lut_progressive_stacks = np.zeros((number_of_stacks, 6))
        lut_stacks_boundaries = np.zeros((number_of_stacks, 4))
        lut_stacks_paths = []

        # sigma nought  and theta initialization
        sigma_nought_equi7_fnames = {}
        theta_equi7_fnames = {}
        # cycle for each stack
        # for stack_key, scene_dict in proc_inputs.stacks_scenes_fnames_orders.items():
        for (stack_progressive_idx, (unique_stack_id, acquisitions_pf_names),) in enumerate(
            input_params_obj.stack_based_processing.stack_composition.items()
        ):

            (
                global_cycle_idx,
                heading_deg,
                rg_swath_idx,
                rg_sub_swath_idx,
                az_swath_idx,
                _,
            ) = decode_unique_acquisition_id_string(unique_stack_id + "_BSL_00")

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
            temp_output_folder_gr = os.path.join(temp_output_folder, "geocoded", unique_stack_id)
            temp_output_folder_e7 = os.path.join(temp_output_folder, "equi7", unique_stack_id)
            os.makedirs(temp_output_folder_gr)
            os.makedirs(temp_output_folder_e7)

            ### load data ( and oversample if requested and if needed )
            try:
                logging.info("AGB: Data loading for stack " + unique_stack_id + "; this may take a while:")

                (data_SLC, master_id, raster_info, raster_info_orig,) = read_and_oversample_data(
                    input_params_obj.dataset_query.L1C_repository,
                    acquisitions_pf_names,
                    conf_params_obj.processing_flags.enable_resampling,
                )

            except Exception as e:
                logging.error("AGB: error during input data reading: " + str(e), exc_info=True)
                raise

            ### load or compute auxiliary data
            try:

                read_ref_h = (
                    not conf_params_obj.processing_flags.apply_calibration_screen
                    and conf_params_obj.processing_flags.DEM_flattening
                )
                read_cal_screens = conf_params_obj.processing_flags.apply_calibration_screen
                geometry_aux_are_present = check_if_geometry_auxiliaries_are_present(
                    input_params_obj.stack_based_processing,
                    unique_stack_id,
                    acquisitions_pf_names,
                    read_ref_h=read_ref_h,
                )

                if conf_params_obj.processing_flags.compute_geometry or not geometry_aux_are_present:

                    # messages for the log:
                    if conf_params_obj.processing_flags.compute_geometry:
                        logging.info("AGB: calling geometry library for stack " + unique_stack_id + "\n")
                        if geometry_aux_are_present:
                            logging.warning("    geometry auxiliaries will be overwritten for stack " + unique_stack_id)
                        else:
                            logging.info("\n")
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
                        input_params_obj.dataset_query.L1C_repository,
                        input_params_obj.stack_based_processing,
                        unique_stack_id,
                        acquisitions_pf_names,
                        master_id,
                        conf_params_obj.processing_flags.enable_resampling,
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
                        input_params_obj.dataset_query.L1C_repository,
                        input_params_obj.stack_based_processing,
                        unique_stack_id,
                        acquisitions_pf_names,
                        master_id,
                        conf_params_obj.processing_flags.enable_resampling,
                        comp_ref_h=read_ref_h,
                        sar_geometry_master=sar_geometry_master,
                    )

                    logging.info("Geometry library: correcting geometry global / local reference")
                    slope = slope - ellipsoid_slope

                    # The kz taking into account the ground slope is generated too
                    kz_ground_slope = kz.copy()

                    for swath_id in kz.keys():
                        kz[swath_id] = (
                            kz[swath_id]
                            * np.sin(off_nadir_angle_rad[master_id])
                            / np.sin(off_nadir_angle_rad[master_id] - ellipsoid_slope)
                        )
                        kz_ground_slope[swath_id] = kz[swath_id] / (
                            1 - np.tan(slope) / np.tan(off_nadir_angle_rad[master_id] - ellipsoid_slope)
                        )
                    off_nadir_angle_rad[master_id] = off_nadir_angle_rad[master_id] - ellipsoid_slope
                    del ellipsoid_slope

                    logging.info("AGB: ...geometry auxiliaries computation done.")

                else:

                    logging.info(
                        "AGB: geometry auxiliaries are provided from user, so they are now being loaded and not computed, for stack "
                        + unique_stack_id
                        + "\n"
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
                        input_params_obj.stack_based_processing,
                        unique_stack_id,
                        acquisitions_pf_names,
                        conf_params_obj.processing_flags.enable_resampling,
                        raster_info_orig,
                        read_ref_h=read_ref_h,
                    )

                    # The kz taking into account the ground slope is generated too
                    # ellipsoid_slope is set to zero
                    kz_ground_slope = kz.copy()
                    for swath_id in kz.keys():
                        kz_ground_slope[swath_id] = kz[swath_id] / (
                            1 - np.tan(slope) / np.tan(off_nadir_angle_rad[master_id])
                        )
                    logging.info("AGB: ...geometry auxiliaries loading done.")

                # read the rest of auxiliaries which are notpart of the geometry library:
                if read_cal_screens:

                    logging.warning("AGB: loading calibration screens \n")

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
                        input_params_obj.stack_based_processing,
                        unique_stack_id,
                        acquisitions_pf_names,
                        conf_params_obj.processing_flags.enable_resampling,
                        raster_info_orig,
                        read_cal_screens=read_cal_screens,
                        read_ecef=False,
                        read_off_nadir=False,
                        read_slope=False,
                        read_kz=False,
                        read_ref_h=False,
                        read_dist=False,
                    )
                    logging.info("...done")

            except Exception as e:
                logging.error(
                    "AGB: error during auxiliary data computation and/or loading: " + str(e), exc_info=True,
                )
                raise

            ### Screen calibration (ground steering)
            try:
                if conf_params_obj.processing_flags.apply_calibration_screen:
                    logging.info("AGB: applying calibration screen...")
                    data_SLC = apply_calibration_screens(
                        data_SLC, raster_info, cal_screens, cal_screens_raster_info, master_id,
                    )
                    logging.info("...done.\n")

                elif conf_params_obj.processing_flags.DEM_flattening:
                    logging.info("AGB: DEM flattening... ")
                    data_SLC = apply_dem_flattening(data_SLC, kz, reference_height, master_id, raster_info)
                    logging.info("...done.\n")

            except Exception as e:
                logging.error(
                    "AGB: error during screen calibration or DEM flattening." + str(e), exc_info=True,
                )
                raise

            ### ground cancellation
            try:

                if space_varying_ground_cancellation:

                    # The map specifying the space-varying height to be emphasized
                    # must be available here. For the moment a map filled with the
                    # value taken from the configuration file is used.
                    current_enhanced_forest_height = slope.copy()
                    current_enhanced_forest_height.fill(conf_params_obj.ground_cancellation.enhanced_forest_height)

                    logging.info("AGB: space-varying ground contribute cancellation...:")
                    SLC_ground_cancelled = ground_cancellation(
                        data_SLC,
                        kz_ground_slope,
                        conf_params_obj.ground_cancellation.multi_master_flag,
                        current_enhanced_forest_height,
                        conf_params_obj.ground_cancellation.equalization_flag,
                        raster_info.resolution_m_slant_rg,
                        off_nadir_angle_rad[master_id],
                        slope,
                    )  # data_SLC: 3 pol (nrg x Naz  x N images); SLC_ground_cancelled: 3 pol (Nrg x Naz x 1 image)
                    del data_SLC, kz, kz_ground_slope, current_enhanced_forest_height
                else:

                    logging.info("AGB: ground contribute cancellation...:")

                    SLC_ground_cancelled = ground_cancellation(
                        data_SLC,
                        kz,
                        conf_params_obj.ground_cancellation.multi_master_flag,
                        conf_params_obj.ground_cancellation.enhanced_forest_height,
                        conf_params_obj.ground_cancellation.equalization_flag,
                        raster_info.resolution_m_slant_rg,
                        off_nadir_angle_rad[master_id],
                        slope,
                    )  # data_SLC: 3 pol (nrg x Naz  x images); SLC_ground_cancelled: 3 pol (Nrg x Naz x 1 image)

                    del data_SLC, kz

            except Exception as e:
                logging.error("AGB: error during ground cancellation: " + str(e), exc_info=True)
                raise
            logging.info("...done.\n")

            ### compute mean look angle
            look_angle_rad = np.nanmean(off_nadir_angle_rad[master_id])
            logging.info("AGB: look angle used is {} [deg] \n".format(np.rad2deg(look_angle_rad)))

            if (
                conf_params_obj.estimate_agb.intermediate_ground_averaging
                > conf_params_obj.estimate_agb.product_resolution / 2
            ):
                sigma_ground_res_m = conf_params_obj.estimate_agb.product_resolution / 2
                logging.warning(
                    '"intermediate ground averaging" cannot be greater than "product resolution/2", setting it to {}'.format(
                        sigma_ground_res_m
                    )
                )
            else:
                sigma_ground_res_m = conf_params_obj.estimate_agb.intermediate_ground_averaging

            if conf_params_obj.processing_flags.multilook_heading_correction:
                sigma_ground_res_m = resolution_heading_correction(sigma_ground_res_m, heading_deg)

            windtm_x = np.int(np.round(sigma_ground_res_m / raster_info.pixel_spacing_az / 2) * 2 + 1)
            windtm_y = np.int(
                np.round(sigma_ground_res_m / (raster_info.pixel_spacing_slant_rg / np.sin(look_angle_rad)) / 2) * 2 + 1
            )

            sub_factor_x = np.int((windtm_x - 1) / 2)
            sub_factor_y = np.int((windtm_y - 1) / 2)

            logging.info("AGB: multilooking of incidence angle...")
            theta_multi_looked_sr = convolve2d(
                off_nadir_angle_rad[master_id] - slope,
                np.ones((windtm_y, windtm_x)) / windtm_y / windtm_x,
                mode="same",
            )
            theta_multi_looked_sr = theta_multi_looked_sr[::sub_factor_y, ::sub_factor_x]
            logging.info("...done.\n")

            sigma_nought_sr = {}
            for pol_name in SLC_ground_cancelled.keys():

                logging.info("AGB: multilooking of ground notched for polarization {}...".format(pol_name))
                SLC_ground_cancelled_multi_look = convolve2d(
                    np.absolute(SLC_ground_cancelled[pol_name]) ** 2,
                    np.ones((windtm_y, windtm_x)) / windtm_y / windtm_x,
                    mode="same",
                )

                # It's not important that data SLC data is a beta or sigma nought: it's important to be proportional to beta or sigma nought.
                if conf_params_obj.processing_flags.input_data_type == "beta_nought":
                    # Data SLC is proportional to beta nought: convert to sigma nought (while subsampling)
                    sigma_nought_sr[pol_name] = radiometric_correction_beta_to_sigma(
                        SLC_ground_cancelled_multi_look[::sub_factor_y, ::sub_factor_x], theta_multi_looked_sr
                    )
                elif conf_params_obj.processing_flags.input_data_type == "sigma_nought":
                    # Data SLC is already proportional to sigma nought: just subsampling
                    sigma_nought_sr[pol_name] = SLC_ground_cancelled_multi_look[::sub_factor_y, ::sub_factor_x]

                else:
                    error_message = "Input data type should be a string called beta_nought or sigma_nought; found {} instead. Check configuration file.".format(
                        conf_params_obj.processing_flags.input_data_type
                    )
                    logging.error(error_message)
                    raise RuntimeError(error_message)

                sigma_nought_sr[pol_name][sigma_nought_sr[pol_name] < 0] = np.NaN

                logging.info("...done.\n")

            del SLC_ground_cancelled_multi_look

            ### saving breakpoints
            if conf_params_obj.processing_flags.save_breakpoints:
                logging.info("AGB: saving breakpoints (in slant range geometry) on " + breakpoints_output_folder)
                post_string = "_SR_" + unique_stack_id

                breakpoint_names = ["ground_cancelled" + post_string]

                save_breakpoints(breakpoints_output_folder, breakpoint_names, [SLC_ground_cancelled])
                logging.info("...done.\n")

            del SLC_ground_cancelled

            az_vec_subs = np.arange(0, raster_info.num_lines, sub_factor_x)
            rg_vec_subs = np.arange(0, raster_info.num_samples, sub_factor_y)

            ### Interpolate it over a regular lat lon grid (with grid data): generate a regular grid for the interpolation by using Max and Min lat lon from the ECEFGRID_LLH (make the vectors a bit longer), for the grid steps use the minimum steps from the ECEFGRID_LLH)
            logging.info(unique_stack_id + ": Geocoding data...")
            try:

                # initialize the geocoding
                min_spacing_m = min(
                    sub_factor_x * raster_info.pixel_spacing_az, sub_factor_y * raster_info.pixel_spacing_slant_rg,
                )
                min_spacing_m = min(min_spacing_m, equi7_sampling_intermediate)

                logging.info("Geocoding spacing set to {} [m]".format(min_spacing_m))

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

                # geocode the sigma nought (three polarizations)
                sigma_nought_gr = {}
                for pol_name in sigma_nought_sr.keys():

                    logging.info("AGB: geocoding the sigma nought for polarization {}...".format(pol_name))
                    sigma_nought_gr[pol_name] = geocoding(
                        sigma_nought_sr[pol_name], lon_in, lat_in, lonMeshed_out, latMeshed_out, valid_values_mask,
                    )

                    logging.info("...done.\n")

                del sigma_nought_sr

                # geocode the theta incidence angle
                logging.info("AGB: Geocoding of incidence angle...")
                theta_multi_looked_gr = geocoding(
                    theta_multi_looked_sr, lon_in, lat_in, lonMeshed_out, latMeshed_out, valid_values_mask,
                )
                logging.info("...done.\n")
                del theta_multi_looked_sr

                logging.info("...done.\n")

            except Exception as e:
                logging.error("AGB: error during geocoding: " + str(e), exc_info=True)
                raise

            ### create GEOTIFF of all the layers (estimation, mask ):
            logging.info(unique_stack_id + ": formatting data to GEOTIFF...")
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

                # geotiff of the sigma nought (three polarizations)
                sigma_nought_ground_fnames = {}
                for pol_name in sigma_nought_gr.keys():
                    sigma_nought_ground_fnames[pol_name] = os.path.join(
                        temp_output_folder_gr, "sigma_nought_" + pol_name + ".tif",
                    )

                    tiff_formatter(
                        sigma_nought_gr[pol_name],
                        sigma_nought_ground_fnames[pol_name],
                        geotransform,
                        gdal_data_format=gdal.GDT_Float32,
                    )

                del sigma_nought_gr

                # geotiff of the theta
                theta_ground_fname = os.path.join(temp_output_folder_gr, "theta.tif")

                tiff_formatter(
                    theta_multi_looked_gr, theta_ground_fname, geotransform, gdal_data_format=gdal.GDT_Float32,
                )

                del theta_multi_looked_gr
                logging.info("...done.\n")

            except Exception as e:
                logging.error("AGB: error during GEOTIFF formatting: " + str(e), exc_info=True)
                raise

            ### formatting data to EQUI7
            logging.info(unique_stack_id + ": formatting into EQUI7 grid...")
            try:

                sigma_nought_equi7_fnames[unique_stack_id] = {}
                # equi7 of the sigma nought (three polarizations)
                for pol_name in sigma_nought_ground_fnames.keys():

                    equi7_sigma_nought_outdir = os.path.join(temp_output_folder_e7, "sigma_nought_" + pol_name)

                    logging.info(
                        "image2equi7grid IN: "
                        + sigma_nought_ground_fnames[pol_name]
                        + " , OUT:"
                        + equi7_sigma_nought_outdir
                    )
                    sigma_nought_equi7_fnames[unique_stack_id][pol_name] = image2equi7grid(
                        e7g_intermediate,
                        sigma_nought_ground_fnames[pol_name],
                        equi7_sigma_nought_outdir,
                        gdal_path=gdal_path,
                        inband=None,
                        subgrid_ids=None,
                        accurate_boundary=False,
                        withtilenamesuffix=False,
                        resampling_type="bilinear",
                        tile_nodata=np.nan,
                    )

                # equi7 of the theta
                equi7_theta_outdir = os.path.join(temp_output_folder_e7, "theta")

                logging.info("image2equi7grid IN: " + theta_ground_fname + " , OUT:" + equi7_theta_outdir)
                theta_equi7_fnames[unique_stack_id] = image2equi7grid(
                    e7g_intermediate,
                    theta_ground_fname,
                    equi7_theta_outdir,
                    gdal_path=gdal_path,
                    inband=None,
                    subgrid_ids=None,
                    accurate_boundary=False,
                    withtilenamesuffix=False,
                    resampling_type="bilinear",
                    tile_nodata=np.nan,
                )

            except Exception as e:
                logging.error("AGB: error during EQUI7 formatting: " + str(e), exc_info=True)
                raise

            if not os.path.exists(theta_equi7_fnames[unique_stack_id][0]):
                error_message = "EQUI7 grid has not been generated, output is absent "
                logging.error(error_message)
                raise RuntimeError(error_message + " :" + theta_equi7_fnames[unique_stack_id][0])

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

            logging.info("...done.\n")

        equi7_initialization = {
            "equi7_sampling_intermediate": equi7_sampling_intermediate,
            "e7g_intermediate": e7g_intermediate,
        }
        # create the CoreProcessingAGB configuration file, starting from the default one
        # Read the default conf file:
        default_coreprocessing_conf_file = os.path.join(
            os.path.dirname(self.configuration_file), "ConfigurationFile_CoreProcessingAGB_Default.xml",
        )

        # update the conf paths:
        # source composition:
        # source[index_obs][index_stack][index_file][path,band_id]
        for index_obs, name in enumerate(conf_params_obj.estimate_agb.residual_function.formula_observables.name):
            if not conf_params_obj.estimate_agb.residual_function.formula_observables.source_paths[index_obs]:

                if name == "neg_sigma0_hh_db":
                    pol_name = "HH"
                    stack_list = []
                    for index_stack, sigma0_pols_dict in enumerate(sigma_nought_equi7_fnames.values()):
                        sigma0_file_names_list = sigma0_pols_dict[pol_name]
                        file_list = []
                        for index_file, sigma0_file_name in enumerate(sigma0_file_names_list):
                            layer_list = [sigma0_file_name, 0]
                            file_list.append(layer_list)
                        stack_list.append(file_list)

                    ## temporary (for tests with tomographic data)
                    # layer_list = [r'C:\Users\macie\Documents\bio_input\aux_for_agb_dev\out_tomo_GGS_50m_RES_100m\slice_30_m\EQUI7_AF050M\E045N048T3\cov_vv_30_m_EQUI7_AF050M_E045N048T3.tif', 0]
                    # file_list = [layer_list]
                    # stack_list = [file_list]

                    conf_params_obj.estimate_agb.residual_function.formula_observables.source_paths[
                        index_obs
                    ] = stack_list
                    conf_params_obj.estimate_agb.residual_function.formula_observables.source_resolution[
                        index_obs
                    ] = sigma_ground_res_m
                    conf_params_obj.estimate_agb.residual_function.formula_observables.source_unit[index_obs] = "none"

                elif name == "neg_sigma0_hv_db":
                    pol_name = "VH"
                    stack_list = []
                    for index_stack, sigma0_pols_dict in enumerate(sigma_nought_equi7_fnames.values()):
                        sigma0_file_names_list = sigma0_pols_dict[pol_name]
                        file_list = []
                        for index_file, sigma0_file_name in enumerate(sigma0_file_names_list):
                            layer_list = [sigma0_file_name, 0]
                            file_list.append(layer_list)
                        stack_list.append(file_list)

                    ## temporary (for tests with tomographic data)
                    # layer_list = [r'C:\Users\macie\Documents\bio_input\aux_for_agb_dev\out_tomo_GGS_50m_RES_100m\slice_30_m\EQUI7_AF050M\E045N048T3\cov_vh_30_m_EQUI7_AF050M_E045N048T3.tif', 0]
                    # file_list = [layer_list]
                    # stack_list = [file_list]

                    conf_params_obj.estimate_agb.residual_function.formula_observables.source_paths[
                        index_obs
                    ] = stack_list
                    conf_params_obj.estimate_agb.residual_function.formula_observables.source_resolution[
                        index_obs
                    ] = sigma_ground_res_m
                    conf_params_obj.estimate_agb.residual_function.formula_observables.source_unit[index_obs] = "none"

                elif name == "neg_sigma0_vv_db":
                    pol_name = "VV"
                    stack_list = []
                    for index_stack, sigma0_pols_dict in enumerate(sigma_nought_equi7_fnames.values()):
                        sigma0_file_names_list = sigma0_pols_dict[pol_name]
                        file_list = []
                        for index_file, sigma0_file_name in enumerate(sigma0_file_names_list):
                            layer_list = [sigma0_file_name, 0]
                            file_list.append(layer_list)
                        stack_list.append(file_list)

                    ## temporary (for tests with tomographic data)
                    # layer_list = [r'C:\Users\macie\Documents\bio_input\aux_for_agb_dev\out_tomo_GGS_50m_RES_100m\slice_30_m\EQUI7_AF050M\E045N048T3\cov_vv_30_m_EQUI7_AF050M_E045N048T3.tif', 0]
                    # file_list = [layer_list]
                    # stack_list = [file_list]

                    conf_params_obj.estimate_agb.residual_function.formula_observables.source_paths[
                        index_obs
                    ] = stack_list
                    conf_params_obj.estimate_agb.residual_function.formula_observables.source_resolution[
                        index_obs
                    ] = sigma_ground_res_m
                    conf_params_obj.estimate_agb.residual_function.formula_observables.source_unit[index_obs] = "none"

                elif (name == "cos_local_db") or (name == "theta_local"):

                    stack_list = []
                    for index_stack, theta_list in enumerate(theta_equi7_fnames.values()):
                        file_list = []
                        for index_file, theta_file_name in enumerate(theta_list):
                            layer_list = [theta_file_name, 0]
                            file_list.append(layer_list)
                        stack_list.append(file_list)

                    ## temporary (for tests with tomographic data)
                    # layer_list = [r'C:\Users\macie\Documents\bio_input\aux_for_agb_dev\out_tomo_GGS_50m_RES_100m\theta\EQUI7_AF050M\E045N048T3\theta_AF050M_E045N048T3.tif', 0]
                    # file_list = [layer_list]
                    # stack_list = [file_list]

                    conf_params_obj.estimate_agb.residual_function.formula_observables.source_paths[
                        index_obs
                    ] = stack_list
                    conf_params_obj.estimate_agb.residual_function.formula_observables.source_resolution[
                        index_obs
                    ] = sigma_ground_res_m
                    conf_params_obj.estimate_agb.residual_function.formula_observables.source_unit[index_obs] = "rad"

                ## temporary (for tests with gedi-like data for calibration)
                # elif name == "agb_1_cal_1km_db":
                #     # layer_list = [r'C:\Users\macie\Documents\bio_input\aux_for_agb_dev\lope_lidar\lidar_agb\EQUI7_AF050M\E045N048T3\lidar_AGB_1km_AF050M_E045N048T3.tif', 0]
                #     # file_list = [layer_list]
                #     # stack_list = [file_list]
                #     stack_list = []
                #     conf_params_default.AGB.residual_function.formula_observables.source_paths[index_obs] = stack_list
                #     conf_params_default.AGB.residual_function.formula_observables.source_resolution[index_obs] = 1000
                #     conf_params_default.AGB.residual_function.formula_observables.source_unit[index_obs] = 't/ha'

                elif name == "agb_1_cal_db":

                    file_list = []
                    for index_file, cal_file_name in enumerate(lut_cal_paths):
                        layer_list = [cal_file_name, 0]
                        file_list.append(layer_list)
                    stack_list = [file_list]

                    conf_params_obj.estimate_agb.residual_function.formula_observables.source_paths[
                        index_obs
                    ] = stack_list
                    conf_params_obj.estimate_agb.residual_function.formula_observables.source_resolution[
                        index_obs
                    ] = 100
                    conf_params_obj.estimate_agb.residual_function.formula_observables.source_unit[index_obs] = "none"

                    ## uncomment if testing with gedi-like data
                    # stack_list = []
                    conf_params_obj.estimate_agb.residual_function.formula_observables.source_paths[
                        index_obs
                    ] = stack_list
                    conf_params_obj.estimate_agb.residual_function.formula_observables.source_resolution[index_obs] = 50
                    conf_params_obj.estimate_agb.residual_function.formula_observables.source_unit[index_obs] = "t/ha"

                ## temporary (for tests with height data)
                # elif (name == "tomo_h_db"):
                #     layer_list = [r'C:\Users\macie\Documents\bio_input\aux_for_agb_dev\out_tomo_fh\BIOMASS_L2_20201118T123101\TOMO_FH\Products\global_FH\EQUI7_AF050M\E045N048T3\FH_EQUI7_AF050M_E045N048T3.tif', 0]
                #     file_list = [layer_list]
                #     stack_list = [file_list]
                #     conf_params_default.AGB.residual_function.formula_observables.source_paths[index_obs] = stack_list
                #     conf_params_default.AGB.residual_function.formula_observables.source_resolution[index_obs] = 100
                #     conf_params_default.AGB.residual_function.formula_observables.source_unit[index_obs] = 'm'

                elif name == "forest_class":

                    file_list = []
                    for index_file, fnf_file_name in enumerate(lut_fnf_paths):
                        layer_list = [fnf_file_name, 0]
                        file_list.append(layer_list)
                    stack_list = [file_list]

                    conf_params_obj.estimate_agb.residual_function.formula_observables.source_paths[
                        index_obs
                    ] = stack_list
                    conf_params_obj.estimate_agb.residual_function.formula_observables.source_resolution[
                        index_obs
                    ] = 100
                    conf_params_obj.estimate_agb.residual_function.formula_observables.source_unit[index_obs] = "none"

        lut_cal = LookupTableAGB(paths=lut_cal_paths, boundaries=lut_cal_boundaries, progressive=None)
        lut_fnf = LookupTableAGB(paths=lut_fnf_paths, boundaries=lut_fnf_boundaries, progressive=None)
        lut_stacks = LookupTableAGB(
            paths=lut_stacks_paths, boundaries=lut_stacks_boundaries, progressive=lut_progressive_stacks,
        )

        # write the input file with the sections needed by the Core Processing AGB APP:
        out_input_file_xml = os.path.join(
            os.path.dirname(input_params_obj.output_specification.output_folder), "Input_File_CoreProcessingAGB.xml"
        )
        input_params_obj.core_processing_agb = fill_core_processing_agb_obj(
            input_params_obj, lut_cal, lut_fnf, lut_stacks
        )
        write_input_file(input_params_obj, out_input_file_xml)

        # write the updatec conf file:
        coreprocessing_configuration_file_xml = os.path.join(
            os.path.dirname(input_params_obj.output_specification.output_folder), "Configuration_File_paths_updated.xml"
        )
        write_configuration_file(conf_params_obj, coreprocessing_configuration_file_xml)

        end_message = "AGB stack-based processing APP ended correctly.\n"
        logging.info(end_message)
        print(end_message)
        ########################## END OF STACK BASED STEPS ######################

        return out_input_file_xml, coreprocessing_configuration_file_xml


class CoreProcessingAGB(Task):
    """CoreProcessingAGB APP

    CoreProcessingAGB is the second of the two sub-APPs called by AboveGroundBiomass processor.

    It computes the AGB product starting from stack inputs generated by "StackBasedProcessingAGB" APP.
    It is automatically called in the default run (see Getting Started  and Tutorial documentation sections),
    or can be manually launched (stand alone) following sequence:
    dataset_query -> StackBasedProcessingAGB -> CoreProcessingAGB

    Attributes
    ----------
    configuration_file : str
        path of the Configuration_File.xml file
        
    Methods
    -------
    run( input_file_path )
        run the CoreProcessingAGB APP
    name : str
        name of the APP

    See Also
    --------
    biopal.dataset_query.dataset_query.dataset_query : it's the first APP to be called in the manual sequence
    StackBasedProcessingAGB : it's the APP which prepares the input for this APP

    Notes
    -----
    Stand alone "CoreProcessingAGB" APP launch. 
    It is possible to launch many istances of "CoreProcessingAGB" APP without need to re launch all the processor.
    
    The app needs a configuration file updated with object paths (updated during "StackBasedProcessingAGB" APP run)
    and can be customized to run different core processings with same stack based pre processing from StackBasedProcessingAGB.

    Examples
    --------
    Manual AGB chain execution

    >>> from biopal.dataset_query.dataset_query import dataset_query
    >>> from biopal.AGB.main_AGB import StackBasedProcessingAGB, CoreProcessingAGB
    >>> dq_obj = dataset_query()
    >>> input_file_up1 = dq_obj.run( input_file )
    >>> sbp_obj = StackBasedProcessingAGB( config_file )
    >>> input_file_up2, config_file_up = sbp_obj.run( input_file_up1 )

    At this point, the CoreProcessingAGB can be launched with different configuration parameters, modifying config_file_up

    >>> agbp_obj = CoreProcessingAGB( config_file_up )
    >>> agbp_obj.run( input_file_up2 )

    - input_file: path of the BioPAL input file
    - input_file_up1: same of input_file with also the "stack_based_processing" section   
    - input_file_up2: same of input_file_up1 with also the "core_processing_agb" section
    - config_file: path of the BioPAL configuration file 
    - config_file_up: same of config_file, with paths updated
    """

    def __init__(self, configuration_file):

        super().__init__(configuration_file)

    def _run(self, input_file):

        # parse the input file and check for the needed sections in it:
        if isinstance(input_file, str):
            check_if_path_exists(input_file, "FILE")
            input_params_obj = parse_input_file(input_file)
        else:
            input_params_obj = input_file
        if (
            input_params_obj.L2_product is None
            or input_params_obj.output_specification is None
            or input_params_obj.dataset_query is None
            or input_params_obj.core_processing_agb is None
        ):
            error_message = [
                "CoreProcessingAGB APP input file should contain at least"
                ' "L2_product", "output_specification", "dataset_query" and "core_processing_agb" sections'
            ]
            raise InvalidInputError(error_message)

        # parse the configuration file and check for the needed sections in it:
        if isinstance(input_file, str):
            check_if_path_exists(self.configuration_file, "FILE")
            conf_params_obj = parse_configuration_file(self.configuration_file)
        else:
            conf_params_obj = self.configuration_file
        if (
            conf_params_obj.processing_flags is None
            or conf_params_obj.ground_cancellation is None
            or conf_params_obj.estimate_agb is None
        ):
            error_message = [
                "Configuration file for  CoreProcessingAGB APP should contain at least"
                ' "processing_flags", "ground_cancellation" and "estimate_agb" sections'
            ]
            raise InvalidInputError(error_message)

        # Initialize the logger, of the APP has been called in a stand-alone call
        if not hasattr(logging.getLoggerClass().root.handlers[0], "baseFilename"):
            start_logging(
                os.path.dirname(input_params_obj.output_specification.output_folder),
                input_params_obj.L2_product,
                "DEBUG",
                app_name="CoreProcessingAGB",
            )
        logging.info("AGB core-processing APP starting\n")

        # get needed parameters from input and configuration files
        gdal_path, _ = set_gdal_paths(conf_params_obj.gdal.gdal_path, conf_params_obj.gdal.gdal_environment_path)

        geographic_boundaries = input_params_obj.stack_based_processing.geographic_boundaries
        lut_cal_paths = input_params_obj.core_processing_agb.lut_cal.paths
        lut_cal = input_params_obj.core_processing_agb.lut_cal.boundaries

        lut_fnf_paths = input_params_obj.core_processing_agb.lut_fnf.paths
        lut_fnf = input_params_obj.core_processing_agb.lut_fnf.boundaries

        lut_stacks_paths = input_params_obj.core_processing_agb.lut_stacks.paths
        lut_stacks_boundaries = input_params_obj.core_processing_agb.lut_stacks.boundaries
        lut_progressive_stacks = input_params_obj.core_processing_agb.lut_stacks.progressive

        ### initialize the equi7 sampling grid
        equi7_sampling_intermediate = choose_equi7_sampling(
            conf_params_obj.estimate_agb.intermediate_ground_averaging,
            conf_params_obj.estimate_agb.intermediate_ground_averaging / 2,
        )
        e7g_intermediate = Equi7Grid(equi7_sampling_intermediate)

        # setting up directories and making sure that preprocessing has been run
        products_folder = os.path.join(input_params_obj.output_specification.output_folder, "Products")
        temp_proc_folder = os.path.join(products_folder, "temp")
        if not (os.path.exists(temp_proc_folder)):
            error_message = '"temp" folder is not present in output: StackBasedProcessingAGB APP should be launched before CoreProcessingAGB '
            logging.error(error_message)
            raise RuntimeError(error_message)

        # preparing folders for final and temporary output
        global_agb_folder = os.path.join(products_folder, "global_AGB")
        temp_agb_folder = os.path.join(temp_proc_folder, "agb_estimation")
        if not os.path.exists(global_agb_folder):
            os.makedirs(global_agb_folder)
        if not os.path.exists(temp_agb_folder):
            os.makedirs(temp_agb_folder)

        ### THIS PART READS THE SETUP

        """
        
        Note: the xml file with setup and the associated reading function must
        be updated so that the formula, observables, and parameters all are defined
        
        
        """
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
            equi7_sampling_intermediate,
            input_params_obj.output_specification.geographic_grid_sampling,
            geographic_boundaries,
            conf_params_obj.estimate_agb,
        )

        #### read formula, parameter and observable defitions
        algorithm_setup = conf_params_obj.estimate_agb
        formula_terms = algorithm_setup.residual_function.formula_terms
        formula_parameters = algorithm_setup.residual_function.formula_parameters
        formula_observables = algorithm_setup.residual_function.formula_observables

        # list with paths for images that will be merged at the end
        parameter_map_pathlists = [[] for formula_parameter_name in formula_parameters.name]

        # convert parameter limits from output units to input units
        # (output limits are easier for the user to define as they are more natural;
        # the code will use transformed limits though)
        for idx, transform in enumerate(formula_parameters.transform):
            if transform != "none":
                formula_parameters.limits[idx][0] = transform_function(
                    formula_parameters.limits[idx][0], [-np.inf, np.inf], transform, do_forward=True,
                )
                formula_parameters.limits[idx][1] = transform_function(
                    formula_parameters.limits[idx][1], [-np.inf, np.inf], transform, do_forward=True,
                )
                formula_parameters.limit_units[idx] = transform + "_" + formula_parameters.limit_units[idx]

        # convert  observable limits to radians if in degrees
        for idx, unit in enumerate(formula_observables.limit_units):
            # if there is just one unit, check if the source unit is rad and limits unit is deg; if so, convet the limits unit
            if unit.lower() == "deg":
                formula_observables.limits[idx][0] = np.deg2rad(formula_observables.limits[idx][0])
                formula_observables.limits[idx][1] = np.deg2rad(formula_observables.limits[idx][1])
                formula_observables.limit_units[idx] = "rad"

        # check source units agains limit units, flag discrepancies:
        for unit_source, unit_limit, name in zip(
            formula_observables.source_unit, formula_observables.limit_units, formula_observables.name,
        ):
            if unit_source.lower() != unit_limit.lower():
                logging.warning(
                    "AGB: source unit and limit unit do not match for parameter {} ({} and {}). Could still be OK. Proceed with caution.".format(
                        name, unit_source, unit_limit
                    )
                )

        # check parameter limit units against associated observable limit units, flag discrepancies
        for unit_limit, associated_observable, name in zip(
            formula_parameters.limit_units, formula_parameters.associated_observable_name, formula_parameters.name,
        ):
            if associated_observable != "none":
                observable_idx = np.where(np.array(formula_observables.name) == associated_observable)[0][0]
                if formula_observables.limit_units[observable_idx].lower() != unit_limit.lower():
                    logging.warning(
                        "AGB: limit units for parameter {} and associated observable {} do not match ({} and {}). However, this could still be OK: the observable source unit is {} and transform is {}. Proceed with caution.".format(
                            name,
                            associated_observable,
                            unit_limit,
                            formula_observables.limit_units[observable_idx],
                            formula_observables.source_unit[observable_idx],
                            formula_observables.transform[observable_idx],
                        )
                    )

        # checking that none of observable and parameter names can be contained within each other
        # (parsing of the formula will fail if this is the case)
        defined_quantity_names = formula_parameters.name + formula_observables.name
        for defined_quantity_name_1 in defined_quantity_names:
            for defined_quantity_name_2 in defined_quantity_names:
                if (defined_quantity_name_1 != defined_quantity_name_2) & (
                    (defined_quantity_name_1.find(defined_quantity_name_2) >= 0)
                    | (defined_quantity_name_2.find(defined_quantity_name_1) >= 0)
                ):
                    logging.error(
                        "AGB: parameter and observable names must be unique and cannot contain each other ({} and {} do not fulfill this requirement).".format(
                            defined_quantity_name_1, defined_quantity_name_2
                        )
                    )

        ### PREPARE EQUI7 GRID
        # output equi7 projection info
        #   here, we assume that the output tile and subtile will be that of the first observable source
        #   that is covered by the current block
        equi7_info_source_path = [
            x[0] for z in formula_observables.source_paths for y in z for x in y if x[0].find("EQUI7") >= 0
        ][0]
        equi7_subtile_name = equi7_info_source_path.split(os.path.sep)[-3:-1][0][6:]
        equi7_tile_name = equi7_info_source_path.split(os.path.sep)[-3:-1][1]

        equi7_subgrid_code = equi7_subtile_name[:2]
        equi7_projection_string = get_projection_from_path(equi7_info_source_path)
        equi7_product = Equi7Grid(geographic_grid_sampling)

        ### PREPARE ADDITIONAL SAMPLING POLYGONS
        # (e.g., GEDI pixels or 3rd party reference data)
        additional_sampling_polygons = []
        # this is temporary; this should be read from an xml-file or something
        test_additional_polygons = False
        if test_additional_polygons:
            # as an example, use 1km resolution cells (e.g., to calibrate with GEDI pixels)
            dd = 1000  # dimension
            east_mesh_flattened, north_mesh_flattened = [
                x.flatten()
                for x in np.meshgrid(
                    np.arange(first_pixel_east, last_pixel_east + 1000, 1000),
                    np.arange(last_pixel_north, first_pixel_north + 1000, 1000),
                )
            ]
            for e, n in zip(east_mesh_flattened, north_mesh_flattened):
                additional_sampling_polygons.append(Polygon([(e, n), (e, n + dd), (e + dd, n + dd), (e + dd, n)]))

        ### PREPARING STACK INFO
        # read acquisition info table
        stack_info_table = lut_progressive_stacks
        stack_info_table_columns = [
            "stack_id",
            "global_cycle_id",
            "heading_degrees",
            "swath_id",
            "subswath_id",
            "azimuth_id",
        ]

        ### PREPARING PARAMETER GRID
        # parameter block mesh and flattened coordinate vectors
        block_mesh_east, block_mesh_north = np.meshgrid(
            np.arange(first_pixel_east, last_pixel_east, block_spacing_east),
            np.arange(first_pixel_north, last_pixel_north, block_spacing_north),
        )
        block_corner_coordinates_east, block_corner_coordinates_north = [
            np.float64(x.flatten()) for x in [block_mesh_east, block_mesh_north]
        ]
        # count the blocks and set up a "finished" flag vector
        number_of_blocks = len(block_corner_coordinates_east.flatten())

        # Compute the block processing order
        block_order = compute_block_processing_order(
            block_corner_coordinates_east,
            block_corner_coordinates_north,
            block_size_east,
            block_size_north,
            lut_cal,  # right now, this uses boundaries but in the future it should be capable of using polygons (including additional_sampling_polygons)
            lut_stacks_boundaries,  # right now, this uses boundaries but in the future it should be capable of using polygons
        )

        ### RUNNING PARAMETER BLOCKS
        # block status vector
        block_status = np.zeros(number_of_blocks)
        # flag explanation
        # -1: critical error
        # 0: block run but skipped for some detected, data-related reason
        # 1: block successfully finished

        for counter_blocks_run, current_block_index in enumerate(block_order):

            # flag indicating whether to skip the current block or not
            # (must be reset upon each iteration)
            skip_current_block = False

            logging.info(
                "AGB: Running block {} out of {} (block ID: {})".format(
                    counter_blocks_run + 1, number_of_blocks, current_block_index
                )
            )

            # %% ### CREATING SAMPLING GRID AND TESTING FOR DATA
            try:

                logging.info("AGB: preparing sampling grid and tabulating data...")

                # extent of current block # [min_east, max_east, min_north, max_north]
                current_block_extents = np.array(
                    [
                        block_corner_coordinates_east[current_block_index],
                        block_corner_coordinates_east[current_block_index] + block_size_east,
                        block_corner_coordinates_north[current_block_index],
                        block_corner_coordinates_north[current_block_index] + block_size_north,
                    ]
                )

                sampling_axis_east = np.arange(current_block_extents[0], current_block_extents[1], sample_spacing_east,)
                sampling_axis_north = np.arange(
                    current_block_extents[2], current_block_extents[3], sample_spacing_north,
                )
                east_mesh_flattened, north_mesh_flattened = [
                    x.flatten() for x in np.meshgrid(sampling_axis_east, sampling_axis_north)
                ]

                # prepare main sampling polygons
                main_sampling_polygons = []
                for e, n in zip(east_mesh_flattened.flatten(), north_mesh_flattened.flatten()):
                    main_sampling_polygons.append(
                        Polygon(
                            [
                                (e, n),
                                (e, n + sample_size_east),
                                (e + sample_size_east, n + sample_size_north),
                                (e + sample_size_north, n),
                            ]
                        )
                    )

                # pixel axes for current block
                pixel_axis_east = np.arange(current_block_extents[0], current_block_extents[1], pixel_size_east)
                pixel_axis_north = np.arange(current_block_extents[2], current_block_extents[3], pixel_size_north)

                # geo transform info
                current_geotransform = [
                    pixel_axis_east[0],
                    pixel_size_east,
                    0,
                    pixel_axis_north[0],
                    0,
                    pixel_size_north,
                ]

                # this should be updated to a more generalised approach based on polygons defining the extent of each image
                (block_has_data, _) = check_block_for_data_and_cal(
                    current_block_extents, lut_stacks_boundaries, lut_cal
                )

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
                    sample_info_table,
                    sample_info_table_columns,
                ) = sample_and_tabulate_data(
                    current_block_extents,
                    pixel_axis_east,
                    pixel_axis_north,
                    main_sampling_polygons + additional_sampling_polygons,
                    block_has_data,
                    stack_info_table,
                    stack_info_table_columns,
                    formula_observables,
                    algorithm_setup.forest_class_observable_name,
                    formula_parameters,
                    number_of_subsets,
                )

                # checking if the tables contain any data
                if observable_table.shape[0] == 0:
                    logging.info("... skipping block #{} due to lack of valid data points".format(current_block_index))
                    skip_current_block = True
                    block_status[counter_blocks_run] = 0

            except Exception as e:
                logging.error(
                    "AGB: error during data sampling and tabulation." + str(e), exc_info=True,
                )
                block_status[counter_blocks_run] = -1
                raise

            if skip_current_block:
                continue

            # %% ### FITTING MODEL TO TABULATED DATA  AND SAVNG

            try:

                logging.info("AGB: fitting model to data...")

                (
                    parameter_tables,
                    space_invariant_parameter_table,
                    space_invariant_parameter_names,
                    space_variant_parameter_table,
                    space_variant_parameter_names,
                ) = fit_formula_to_random_subsets(
                    formula_terms,
                    number_of_subsets,
                    observable_table,
                    observable_names,
                    identifier_table,
                    identifier_names,
                    parameter_position_table,
                    formula_parameters.name,
                    parameter_tables,
                    parameter_table_columns,
                    formula_parameters.parameter_variabilities,
                    algorithm_setup.fraction_of_cal_per_test,
                    algorithm_setup.fraction_of_roi_per_test,
                    algorithm_setup.min_number_of_cals_per_test,
                    algorithm_setup.min_number_of_rois_per_test,
                    algorithm_setup.transfer_function_name,
                )

                # checking if the tables contain any data
                if len(space_invariant_parameter_table) == 0:
                    logging.info("... skipping block #{} due to no subsets found.".format(current_block_index))
                    skip_current_block = True
                    block_status[counter_blocks_run] = 0
                else:
                    logging.info("AGB: saving data to tables.")
                    # save data
                    parameter_property_names = parameter_table_columns[0][-number_of_subsets - 4 :]
                    line_number_string = ["row"]

                    # select formatting for the output tables
                    table_delimiter = "\t"
                    table_precision = 3
                    table_column_width = 25
                    data_type_lut = [
                        [
                            line_number_string,
                            sample_info_table_columns,
                            identifier_names,
                            parameter_position_table_columns,
                            parameter_property_names,
                            observable_names,
                            formula_parameters.name,
                        ],
                        ["d", "d", "d", "d", "f", "f", "f"],
                    ]

                    # save results table (with parameter estimates replicated across samples etc, if necessary; no indices connecting these to parameter tables)
                    curr_path = os.path.join(temp_agb_folder, "results_table_block_{}.txt".format(current_block_index),)
                    curr_table = np.column_stack(
                        (
                            np.arange(observable_table.shape[0]),
                            sample_info_table,
                            identifier_table,
                            observable_table,
                            space_invariant_parameter_table,
                            space_variant_parameter_table,
                            parameter_position_table,
                        )
                    )
                    curr_column_names = (
                        line_number_string
                        + sample_info_table_columns
                        + identifier_names
                        + observable_names
                        + space_invariant_parameter_names
                        + space_variant_parameter_names
                        + parameter_position_table_columns
                    )
                    save_human_readable_table(
                        curr_path,
                        curr_table,
                        curr_column_names,
                        data_type_lut,
                        table_delimiter,
                        table_precision,
                        table_column_width,
                    )
                    curr_dir = os.path.join(temp_agb_folder, "parameter_estimates_subsets")
                    if not os.path.exists(curr_dir):
                        os.makedirs(curr_dir)
                    # save parameter tables (minimalistic tables with only the most necessary info)
                    for parameter_idx, parameter_name in enumerate(formula_parameters.name):
                        curr_path = os.path.join(
                            curr_dir, "parameter_{}_table_block_{}.txt".format(parameter_name, current_block_index),
                        )
                        curr_table = np.column_stack(
                            (np.arange(parameter_tables[parameter_idx].shape[0]), parameter_tables[parameter_idx],)
                        )
                        curr_column_names = np.concatenate(
                            (np.array(line_number_string), parameter_table_columns[parameter_idx],)
                        )
                        save_human_readable_table(
                            curr_path,
                            curr_table,
                            curr_column_names,
                            data_type_lut,
                            table_delimiter,
                            table_precision,
                            table_column_width,
                        )

            except Exception as e:
                logging.error(
                    "AGB: error during parameter estimation or data saving." + str(e), exc_info=True,
                )
                block_status[counter_blocks_run] = -1
                raise

            if skip_current_block:
                continue

            # %% ### READING IMAGE DATA
            try:

                logging.info("AGB: reading data images.")

                # this needs to be updated in line with the tabulation function
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
                    observable_names,
                    formula_observables.source_paths,
                    formula_observables.transform,
                    formula_observables.averaging_method,
                    formula_observables.limits,
                    formula_observables.is_required,
                    [[x, 0] for x in lut_fnf_paths],
                    lut_fnf,
                    identifier_table[:, 2],
                    identifier_table[:, 1],
                    space_invariant_parameter_table,
                    space_invariant_parameter_names,
                    mask_out_area_outside_block=True,
                )

                # skipping to next if all parameters and/or forest class data are missing
                if np.all(forest_class_3d == 0) | np.all(
                    np.array(
                        [
                            np.all(np.isnan(curr_space_invariant_parameter_3d))
                            for curr_space_invariant_parameter_3d in space_invariant_parameters_3d
                        ]
                    )
                ):
                    logging.info("... skipping block #{} due to no image data read.".format(current_block_index))
                    skip_current_block = True
                    block_status[counter_blocks_run] = 0

            except Exception as e:
                logging.error("AGB: error during image reading." + str(e), exc_info=True)
                block_status[counter_blocks_run] = -1
                raise

            if skip_current_block:
                continue

            # %% ### ESTIMATING SPACE-VARIANT PARAMTERS

            try:

                logging.info("AGB: creating space variant parameter images.")

                # detect observables without data (used to remove formula terms that would generate nans)
                observables_with_partial_data = [
                    observable_name
                    for is_nan, observable_name in zip(np.any(np.isnan(observable_table), axis=0), observable_names)
                    if is_nan
                ]
                terms_with_nan_observables = np.any(
                    match_string_lists(formula_terms.string, observables_with_partial_data) >= 0, axis=1,
                )

                if np.any(terms_with_nan_observables):
                    logging.warning(
                        "AGB: skipping formula terms: {} due to lack of useful data for observables: {}.".format(
                            ", ".join(["%d" % (ii + 1) for ii in np.where(terms_with_nan_observables)[0]]),
                            ", ".join(observables_with_partial_data),
                        )
                    )

                #### select relevant formula and weights for this step
                terms_with_zero_weight = np.array(formula_terms.formula_weights.step2) == 0
                if np.any(terms_with_zero_weight):
                    logging.warning(
                        "AGB: skipping formula terms: {} due to zero weights.".format(
                            ", ".join(["%d" % (ii + 1) for ii in np.where(terms_with_zero_weight)[0]])
                        )
                    )

                terms_to_take = ~(terms_with_nan_observables | terms_with_zero_weight)
                formula = subset_iterable(formula_terms.string, terms_to_take)
                formula_weights = subset_iterable(formula_terms.formula_weights.step2, terms_to_take)

                # observables_with_partial_data = [observable_name for is_nan,observable_name in zip(np.any(np.isnan(observable_table),axis=0),observable_names) if is_nan]

                # take out the observables that are in formula and not among space variant parameters
                parameters_for_mapping = np.any(
                    match_string_lists(formula, formula_parameters.name) >= 0, axis=0
                ) & ~np.any(match_string_lists(space_invariant_parameter_names, formula_parameters.name) >= 0, axis=0,)

                observables_for_mapping = np.all(
                    match_string_lists(observable_names, observables_with_partial_data) == -1, axis=1,
                )

                if np.sum(parameters_for_mapping) > 1:
                    logging.error("AGB: the current implementation requires only one space-variant parameter.")

                (space_variant_parameters_3d, space_variant_parameters_3d_names,) = map_space_variant_parameters(
                    formula,
                    formula_weights,
                    forest_class_3d,
                    subset_iterable(observables_3d, observables_for_mapping),
                    subset_iterable(observables_3d_names, observables_for_mapping),
                    space_invariant_parameters_3d,
                    space_invariant_parameters_3d_names,
                    identifiers_3d,
                    identifiers_3d_names,
                    subset_iterable(formula_parameters.name, parameters_for_mapping, False),
                    subset_iterable(formula_parameters.parameter_variabilities, parameters_for_mapping, False,),
                    subset_iterable(formula_parameters.limits, parameters_for_mapping, False),
                    algorithm_setup.transfer_function_name,
                )

                # skipping to next if all parameters and/or forest class data are missing
                if np.all(
                    np.array(
                        [
                            np.all(np.isnan(curr_space_variant_parameter_3d))
                            for curr_space_variant_parameter_3d in space_variant_parameters_3d
                        ]
                    )
                ):
                    logging.info(
                        "... skipping block #{} due to no space variant parameter maps created".format(
                            current_block_index
                        )
                    )
                    skip_current_block = True
                    block_status[counter_blocks_run] = 0

            except Exception as e:
                logging.error(
                    "AGB: error during space variant parameter mapping." + str(e), exc_info=True,
                )
                block_status[counter_blocks_run] = -1
                raise

            if skip_current_block:
                continue

            # %% ### ESTIMATING OTHER PARAMETERS

            try:

                # this is where other things such as error parameters are estimated
                logging.info("AGB: estimating other parameters to be implemented... (includes estimation of error)")

                # for now, the three error elements are hard coded
                additional_parameters_3d = []
                additional_parameters_3d_names = []

            except Exception as e:
                logging.error(
                    "AGB: error during estimation of other parameters." + str(e), exc_info=True,
                )
                block_status[counter_blocks_run] = -1
                raise

            if skip_current_block:
                continue

            # %% ### SAVING SPACE-VARIANT PARAMTERS:

            try:
                logging.info("AGB: saving selected parameters to images.")

                # combine all available mapped parameters
                all_parameters_3d_names = (
                    space_invariant_parameters_3d_names
                    + space_variant_parameters_3d_names
                    + additional_parameters_3d_names
                )
                all_parameters_3d = (
                    space_invariant_parameters_3d + space_variant_parameters_3d + additional_parameters_3d
                )

                for parameter_idx, parameter_name in enumerate(formula_parameters.name):
                    if formula_parameters.save_as_map[parameter_idx]:

                        current_image_to_write = all_parameters_3d[
                            np.where(np.array(all_parameters_3d_names) == parameter_name)[0][0]
                        ]

                        current_file_path = os.path.join(temp_agb_folder, parameter_name)

                        if formula_parameters.transform[parameter_idx] != "none":
                            current_image_to_write = transform_function(
                                current_image_to_write,
                                formula_parameters.limits[parameter_idx],
                                formula_parameters.transform[parameter_idx],
                                do_forward=False,
                            )
                            current_file_path += "_backtransf"

                        current_file_path += (
                            "_block_%d" % (current_block_index) + ".tif"
                        )  # + '_'+equi7_subtile_name+'_'+equi7_tile_name

                        # save tiff file
                        current_file_path = tiff_formatter(
                            [x for x in current_image_to_write.transpose([2, 0, 1])],
                            current_file_path,
                            current_geotransform,
                            gdal_data_format=gdal.GDT_Float32,
                            projection=equi7_projection_string,
                            multi_layers_tiff=True,
                        )

                        # get extents of the current EQUI7 tile
                        # [(left, lower), (right, upper)]
                        try:
                            lon_min, lat_min = getattr(e7g_intermediate, equi7_subgrid_code).xy2lonlat(
                                min(pixel_axis_east), min(pixel_axis_north)
                            )
                            lon_max, lat_max = getattr(e7g_intermediate, equi7_subgrid_code).xy2lonlat(
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
                            gdal_path=gdal_path,
                            ftiles=equi7_product.search_tiles_in_roi(bbox=[(lon_min, lat_min), (lon_max, lat_max)]),
                            accurate_boundary=False,
                            withtilenamesuffix=False,
                            tile_nodata=np.nan,
                        )

                        if formula_parameters.associated_observable_name[parameter_idx] != "none":
                            current_position_in_observable_vector = np.where(
                                match_string_lists(
                                    observable_names, [formula_parameters.associated_observable_name[parameter_idx]],
                                ).flatten()
                                >= 0
                            )[0][0]
                            # current_resolution = 50 # placeholder
                            # current_unit = formula_parameters.units[parameter_idx]
                            if len(formula_observables.source_paths[current_position_in_observable_vector]) == 0:
                                formula_observables.source_paths[current_position_in_observable_vector].append(
                                    [[output_equi7_file_path[0], 0]]
                                )
                            elif len(formula_observables.source_paths[current_position_in_observable_vector]) == 1:
                                formula_observables.source_paths[current_position_in_observable_vector][0].append(
                                    [output_equi7_file_path[0], 0]
                                )
                        parameter_map_pathlists[parameter_idx].append([output_equi7_file_path[0], 0])

                        # lut_cal_paths.append(output_file_path)
                lut_cal = np.row_stack(
                    (lut_cal, np.concatenate((current_block_extents[np.array([0, 1, 3, 2])], np.zeros(1))),)
                )

                skip_current_block = False
                block_status[counter_blocks_run] = 1

            except Exception as e:
                logging.error("AGB: error during saving of maps." + str(e), exc_info=True)
                raise

        # %% FINAL MERGING OF THE IMAGES
        try:
            logging.info("AGB: creating wall-to-wall maps...")
            for parameter_idx, parameter_name in enumerate(formula_parameters.name):
                if formula_parameters.save_as_map[parameter_idx]:

                    tiles_to_save = {}
                    for current_source in parameter_map_pathlists[parameter_idx]:
                        subtile_name = current_source[0].split(os.path.sep)[-3:-1][0][6:]
                        tile_name = current_source[0].split(os.path.sep)[-3:-1][1]

                        if tile_name in tiles_to_save.keys():
                            tiles_to_save[tile_name].append(current_source)
                        else:
                            tiles_to_save[tile_name] = [current_source]

                    for tile_name, tile_sources in tiles_to_save.items():

                        driver = gdal.Open(tile_sources[0][0], GA_ReadOnly)
                        data_merged = np.nan * np.zeros((driver.RasterYSize, driver.RasterXSize))
                        driver = None

                        for idx_tile, current_tile_source in enumerate(tile_sources):

                            # merging current tile blocks togeter:
                            driver = gdal.Open(current_tile_source[0], GA_ReadOnly)

                            if idx_tile == 0:
                                geotransform_out = driver.GetGeoTransform()
                            elif geotransform_out != driver.GetGeoTransform():
                                err_str = "Same equi7 tiles cannot have different geotrasform"
                                logging.error(err_str)
                                raise ValueError(err_str)

                            data_merged = np.nanmean(
                                np.dstack(
                                    (data_merged, driver.GetRasterBand(current_tile_source[1] + 1).ReadAsArray(),)
                                ),
                                axis=2,
                            )

                            driver = None

                        current_merged_file_path = os.path.join(
                            global_agb_folder, equi7_subtile_name, tile_name, parameter_name,
                        )

                        if formula_parameters.transform[parameter_idx] != "none":
                            # current_position_in_observable_vector = parameter_position_in_observable_vector[parameter_idx]
                            # if current_position_in_observable_vector != -1:
                            current_merged_file_path += "_backtransf_"
                        current_merged_file_path += ".tif"

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
                            "... successful for parameter '{}' and filename '{}'".format(
                                parameter_name, current_merged_file_path
                            )
                        )

            end_message = "AGB core-processing APP ended correctly.\n"
            logging.info(end_message)
            print(end_message)

        except Exception as e:
            logging.error(
                "AGB: core-processing APP error during creation of wall-to-wall maps." + str(e), exc_info=True,
            )
            raise


def fill_core_processing_agb_obj(input_params_obj, lut_cal, lut_fnf, lut_stacks):

    """
    Internal function called by the StackBasedProcessingAGB APP:
        
        the StackBasedProcessingAGB APP fills the structure
        to be written into the xml input file for the next APP, which is the 
        core processing for the AGB.
        The returned object "core_processing_agb_obj" contains three lookup tables with paths, boundaries 
        and other specific parameters, of the products (ground notched stacks, 
        forest mask, calibrations) computed by the APP itself.
        
        Usage of the returned object:
            The returned object can be added to the input_params_obj and it 
            can be written to disk if needed, as (core_processing_agb_obj is 
                                                  overwritten by this command, if already present):
                - input_params_obj.core_processing_agb_obj = core_processing_agb_obj
                - write_input_file(input_params_obj, input_file_xml)
    """

    core_processing_agb_obj = core_processing_agb(lut_cal, lut_fnf, lut_stacks,)

    return core_processing_agb_obj
