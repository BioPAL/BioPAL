# SPDX-FileCopyrightText: BioPAL <biopal@esa.int>
# SPDX-License-Identifier: MIT

import os
import numpy as np
import logging
import shutil
from gdalconst import GA_ReadOnly
from osgeo import gdal
from equi7grid.equi7grid import Equi7Grid
from equi7grid.image2equi7grid import image2equi7grid
from copy import deepcopy
from arepytools.timing.precisedatetime import PreciseDateTime

# biomassL2 processor imports
from biopal.fd.processing_FD import alg_wishart_SU
from biopal.data_operations.data_operations import (
    read_and_oversample_data,
    read_and_oversample_aux_data,
    fnf_equi7_load_filter_equi7format,
    fnf_tandemx_load_filter_equi7format,
    apply_dem_flattening,
    get_equi7_tiff_names,
)
from biopal.utility.utility_functions import (
    Task,
    choose_equi7_sampling,
    check_if_path_exists,
    check_if_geometry_auxiliaries_are_present,
    check_fnf_folder_format,
    get_min_time_stamp_repository,
    resolution_heading_correction,
    decode_unique_acquisition_id_string,
    save_breakpoints,
    set_gdal_paths,
    get_data_time_stamp,
    check_equi7_mask_coverage,
)
from biopal.statistics.utility_statistics import (
    main_covariance_estimation_SR,
    covariance_matrix_mat2vec,
    covariance_matrix_vec2mat,
)

from biopal.geocoding.geocoding import (
    geocoding,
    geocoding_init,
)
from biopal.io.xml_io import (
    parse_input_file,
    parse_configuration_file,
)

from biopal.io.data_io import tiff_formatter
from biopal.screen_calibration.screen_calibration import apply_calibration_screens
from biopal.geometry.utility_geometry import compute_and_oversample_geometry_auxiliaries
from biopal.ground_cancellation.ground_cancellation import ground_cancellation


class ForestDisturbance(Task):
    """
    FD main APP "ForestDisturbance" (see BioPAL README.md to launch) is composed by 
    two sub APPS: 
        
    the first APP "initializeAlgorithmFD"  organizes the input dataSet by 
    grouping all the different temporal global cycles of each stack with same 
    nominal geometry
    
    the second APP "CoreProcessingFD" is called with a for-loop for each 
    nominal geometry (each stack) and performs the FD detection algorithm    
    """

    def __init__(
        self, configuration_file,
    ):
        super().__init__(configuration_file)

    def _run(self, input_file):

        # 1) Prepare initialization APP
        # organize the input dataSet by grouping all the temporal global cycles of each stack
        initialize_fd_obj = initializeAlgorithmFD(self.configuration_file)

        # 1) Run the initialization APP
        (
            cycles_composition,  # dictionary
            time_tag_mjd_initial,
            equi7_sampling,
            proc_inputs,
            proc_conf,
        ) = initialize_fd_obj.run(input_file)

        # 2) Main APP: CoreProcessingFD
        # cycle over each nominal geometry stack (containing all its global cycles)

        final_forest_mask_data = {}  # one mask for each equi7 tile, it is updated at each stack-cycle
        for stack_idx, (nominal_geometry_stack_id, global_cycle_dict) in enumerate(cycles_composition.items()):

            # 2) Prepare Main APP: CoreProcessingFD
            core_processing_obj = CoreProcessingFD(
                self.configuration_file,
                cycles_composition,
                stack_idx,
                nominal_geometry_stack_id,
                global_cycle_dict,
                time_tag_mjd_initial,
                final_forest_mask_data,  # optional: this comes from previous iteration
                equi7_sampling,  # optional
                proc_inputs,  # optional
                proc_conf,  # optional
            )

            # 2) Run Main APP: CoreProcessingFD
            (
                final_forest_mask_data,  # this goes in input to next iteration
                temp_output_folder,
            ) = core_processing_obj.run(input_file)

        # Ending the ForestDisturbance APP
        if proc_conf.processing_flags.delete_temporary_files:
            try:
                shutil.rmtree(temp_output_folder)
            except:
                pass
        logging.info("FD processor ended correctly.\n")


class initializeAlgorithmFD(Task):
    """
    Organizes the input dataSet by grouping all the different temporal global 
    cycles of each stack with same nominal geometry, preparing the inputs for 
    the "CoreProcessingFD" APP
    """

    def __init__(self, configuration_file):

        super().__init__(configuration_file)

    def _run(self, input_file):

        ########################## INITIAL STEPS ##############################

        logging.info("FD: Reading chains configuration files")
        check_if_path_exists(self.configuration_file, "FILE")
        proc_conf = parse_configuration_file(self.configuration_file)
        proc_inputs = parse_input_file(input_file)

        ### managing output folders:
        products_folder = os.path.join(proc_inputs.output_specification.output_folder, "Products")
        if proc_conf.processing_flags.save_breakpoints:
            breakpoints_output_folder = os.path.join(products_folder, "breakpoints")
            logging.info("FD: Breakpoints will be saved into: " + breakpoints_output_folder)
            os.makedirs(breakpoints_output_folder)

        temp_output_folder = os.path.join(products_folder, "temp")
        logging.info("FD: Temporary data folder:" + temp_output_folder + "\n")

        os.makedirs(temp_output_folder)

        ### grouping temoral global cycles for each stack
        cycles_composition = (
            {}
        )  # dict which groupes togheter all the differente global cycles for stacks with same nominal geometry
        # cycles_composition[nominal_geometry_stack_id][global_cycle_number]
        logging.info('FD: grouping togheter the stacks which have same "nominal geometry" in different global cycles:')
        for idx, (unique_stack_id, unique_acq_pf_names) in enumerate(
            proc_inputs.stack_based_processing.stack_composition.items()
        ):
            nominal_geometry_stack_id = unique_stack_id[
                6:
            ]  # remove the global cycle indication, maintain the geometry indications
            global_cycle_idx = int(unique_stack_id[3:5])  # \ the cycle number of curent uniqie stack

            if not nominal_geometry_stack_id in cycles_composition.keys():
                cycles_composition[nominal_geometry_stack_id] = {}

            cycles_composition[nominal_geometry_stack_id][global_cycle_idx] = unique_acq_pf_names
        logging.info(
            "    for the #{} stacks in input, #{} different nominal geometries have been found \n".format(
                len(proc_inputs.stack_based_processing.stack_composition), len(cycles_composition)
            )
        )

        ### initialize the equi7 sampling grid
        equi7_sampling = choose_equi7_sampling(
            proc_conf.change_detection_fd.product_resolution, proc_inputs.output_specification.geographic_grid_sampling
        )
        logging.info("EQUI7 Grid sampling used: {}".format(equi7_sampling))

        ### get temporal date time of the input data (get the minimum date from all the stacks)
        time_tag_mjd_initial = get_min_time_stamp_repository(
            proc_inputs.dataset_query.L1C_repository, proc_inputs.stack_based_processing.stack_composition
        )

        ########################## INITIAL STEPS END #############################

        return cycles_composition, time_tag_mjd_initial, equi7_sampling, proc_inputs, proc_conf


class CoreProcessingFD(Task):
    """
    Performs stack-based disturbance algorithm cycling over
    all the different temporal global cycles in the input stack.
    It is called externally from the ForestDisturbance Main APP, once for each
    nominal geometry (stack)
    """

    def __init__(
        self,
        configuration_file,
        cycles_composition,
        stack_idx,
        nominal_geometry_stack_id,
        global_cycle_dict,
        time_tag_mjd_initial,
        final_forest_mask_data=None,
        equi7_sampling=None,
        proc_inputs=None,
        proc_conf=None,
    ):

        super().__init__(configuration_file)

        self.cycles_composition = cycles_composition
        self.stack_idx = stack_idx
        self.nominal_geometry_stack_id = nominal_geometry_stack_id
        self.global_cycle_dict = global_cycle_dict
        self.time_tag_mjd_initial = time_tag_mjd_initial
        self.final_forest_mask_data = final_forest_mask_data
        self.equi7_sampling = equi7_sampling
        self.proc_inputs = proc_inputs
        self.proc_conf = proc_conf

    def check_auxiliaries(self, input_file):
        if self.final_forest_mask_data is None:
            self.final_forest_mask_data = {}
        if self.proc_inputs is None:
            self.proc_inputs = parse_input_file(input_file)
        if self.proc_conf is None:
            self.proc_conf = parse_configuration_file(self.configuration_file)
        if self.equi7_sampling is None:
            self.equi7_sampling = choose_equi7_sampling(
                self.proc_conf.change_detection_fd.product_resolution,
                self.proc_inputs.output_specification.geographic_grid_sampling,
            )

    def _run(self, input_file):

        self.check_auxiliaries(input_file)

        # get needed parameters from input and configuration files
        geographic_boundaries = self.proc_inputs.stack_based_processing.geographic_boundaries

        gdal_path, _ = set_gdal_paths(self.proc_conf.gdal.gdal_path, self.proc_conf.gdal.gdal_environment_path)

        logging.info("FD stack-based processing APP starting\n")

        ########################## INITIALIZATIONS ############################
        number_of_pols = 3

        products_folder = os.path.join(self.proc_inputs.output_specification.output_folder, "Products")
        temp_output_folder = os.path.join(products_folder, "temp")
        if not os.path.exists(temp_output_folder):
            os.makedirs(temp_output_folder)
        if self.proc_conf.processing_flags.save_breakpoints:
            breakpoints_output_folder = os.path.join(products_folder, "breakpoints")
            if not os.path.exists(breakpoints_output_folder):
                os.makedirs(breakpoints_output_folder)

        ### initialize the equi7 sampling grid
        e7g = Equi7Grid(self.equi7_sampling)

        _, heading_deg, rg_swath_idx, rg_sub_swath_idx, az_swath_idx, _ = decode_unique_acquisition_id_string(
            next(iter(self.global_cycle_dict.values()))[0]
        )
        logging.info("FD: computing disturbance for following nominal geometry " + self.nominal_geometry_stack_id + ":")
        logging.info(
            "Heading = {}, RG swath = {} RG sub-swath = {}, AZ swath = {}".format(
                heading_deg, rg_swath_idx, rg_sub_swath_idx, az_swath_idx
            )
        )

        number_of_global_cycles = len(self.global_cycle_dict.keys())
        logging.info("FD: #{} global cycles found for the current nominal geometry".format(number_of_global_cycles))
        ######################### INITIALIZATIONS END #########################

        # cycle over all global cycles of current nominal geometry (current stack)
        for time_step_idx, (global_cycle_idx, uniqie_acq_ids_all_cycles_list) in enumerate(
            self.global_cycle_dict.items()
        ):

            time_tag_mjd_curr = get_data_time_stamp(
                self.proc_inputs.dataset_query.L1C_repository, uniqie_acq_ids_all_cycles_list[0]
            )

            # reconstruct current unique_stack_id
            global_cycle_idx_str = str(global_cycle_idx)
            if len(global_cycle_idx_str) == 1:
                global_cycle_idx_str = "0" + global_cycle_idx_str
            unique_stack_id = "GC_" + global_cycle_idx_str + "_" + self.nominal_geometry_stack_id

            if time_step_idx == 0:
                global_cycle_idx_prev_str = str(global_cycle_idx - 1)
            else:
                global_cycle_idx_prev_str = str([*self.global_cycle_dict.keys()][time_step_idx - 1])
                logging.debug(
                    "time_step_idx {}, global_cycle_idx_prev_str {} ".format(time_step_idx, global_cycle_idx_prev_str)
                )

            if len(global_cycle_idx_prev_str) == 1:
                global_cycle_idx_prev_str = "0" + global_cycle_idx_prev_str

            unique_stack_id_prev = "GC_" + global_cycle_idx_prev_str + "_" + self.nominal_geometry_stack_id

            logging.info(
                "FD: disturbance computation step #{} of #{}; current global cycle is ".format(
                    time_step_idx + 1, number_of_global_cycles
                )
                + '"GC_'
                + global_cycle_idx_str
                + '"'
                + " current geometry is "
                + self.nominal_geometry_stack_id
            )

            # make temporary sub-directories
            temp_output_folder_sr = os.path.join(temp_output_folder, "slantRange", unique_stack_id)
            temp_output_folder_gr = os.path.join(temp_output_folder, "geocoded", unique_stack_id)
            temp_output_folder_e7 = os.path.join(temp_output_folder, "equi7", unique_stack_id)
            os.makedirs(temp_output_folder_sr)
            os.makedirs(temp_output_folder_gr)
            os.makedirs(temp_output_folder_e7)

            ### load data ( and oversample if requested and if needed )
            try:
                logging.info("FD: input data reading:...")

                beta0_calibrated, master_id, raster_info, raster_info_orig = read_and_oversample_data(
                    self.proc_inputs.dataset_query.L1C_repository,
                    uniqie_acq_ids_all_cycles_list,
                    self.proc_conf.processing_flags.enable_resampling,
                )

            except Exception as e:
                logging.error("FD: error during input data reading: " + str(e), exc_info=True)
                raise
            logging.info("...done.\n")

            ### load or compute auxiliary data
            try:

                read_ref_h = (
                    not self.proc_conf.processing_flags.apply_calibration_screen
                    and self.proc_conf.processing_flags.DEM_flattening
                )
                read_cal_screens = self.proc_conf.processing_flags.apply_calibration_screen
                geometry_aux_are_present = check_if_geometry_auxiliaries_are_present(
                    self.proc_inputs.stack_based_processing,
                    unique_stack_id,
                    uniqie_acq_ids_all_cycles_list,
                    read_ref_h=read_ref_h,
                )

                if self.proc_conf.processing_flags.compute_geometry or not geometry_aux_are_present:

                    # messages for the log:
                    if self.proc_conf.processing_flags.compute_geometry:
                        logging.info("FD: calling geometry library for stack " + unique_stack_id + "\n")
                        if geometry_aux_are_present:
                            logging.warning("    geometry auxiliaries will be overwritten for stack " + unique_stack_id)
                        else:
                            logging.info("\n")
                    else:
                        logging.warning(
                            'FD: calling geometry library since AuxiliaryProductsFolder "Geometry" is empty or not complete \n'
                        )

                    _, _, ellipsoid_slope, _, _, _, sar_geometry_master = compute_and_oversample_geometry_auxiliaries(
                        self.proc_inputs.dataset_query.L1C_repository,
                        self.proc_inputs.stack_based_processing,
                        unique_stack_id,
                        uniqie_acq_ids_all_cycles_list,
                        master_id,
                        self.proc_conf.processing_flags.enable_resampling,
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
                        self.proc_inputs.dataset_query.L1C_repository,
                        self.proc_inputs.stack_based_processing,
                        unique_stack_id,
                        uniqie_acq_ids_all_cycles_list,
                        master_id,
                        self.proc_conf.processing_flags.enable_resampling,
                        sar_geometry_master=sar_geometry_master,
                    )

                    logging.info("Geometry library: correcting geometry global / local reference")
                    slope = slope - ellipsoid_slope
                    for swath_id in kz.keys():
                        kz[swath_id] = (
                            kz[swath_id]
                            * np.sin(off_nadir_angle_rad[master_id])
                            / np.sin(off_nadir_angle_rad[master_id] - ellipsoid_slope)
                        )
                    off_nadir_angle_rad[master_id] = off_nadir_angle_rad[master_id] - ellipsoid_slope
                    del ellipsoid_slope

                    logging.info("FD: ...geometry auxiliaries computation done.")

                else:

                    logging.info(
                        "FD: geometry auxiliaries are provided from user, so they are now being loaded and not computed, for stack "
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
                        self.proc_inputs.stack_based_processing,
                        unique_stack_id,
                        uniqie_acq_ids_all_cycles_list,
                        self.proc_conf.processing_flags.enable_resampling,
                        raster_info_orig,
                    )
                    logging.info("FD: ...geometry auxiliaries loading done.")

                # read the rest of auxiliaries which are notpart of the geometry library:
                if read_cal_screens:

                    logging.warning("FD: loading calibration screens \n")

                    _, _, _, _, _, _, _, cal_screens, cal_screens_raster_info, _, _, _ = read_and_oversample_aux_data(
                        self.proc_inputs.stack_based_processing,
                        unique_stack_id,
                        uniqie_acq_ids_all_cycles_list,
                        self.proc_conf.processing_flags.enable_resampling,
                        raster_info_orig,
                        read_cal_screens=read_cal_screens,
                        read_ecef=False,
                        read_off_nadir=False,
                        read_slope=False,
                        read_kz=False,
                        read_ref_h=False,
                        read_dist=False,
                    )
                    logging.info("...done.\n")

            except Exception as e:
                logging.error("FD: error during auxiliary data computation and/or loading: " + str(e), exc_info=True)
                raise

            ### Screen calibration (ground steering)
            try:
                if self.proc_conf.processing_flags.apply_calibration_screen:
                    logging.info("FH: applying calibration screen...")
                    beta0_calibrated = apply_calibration_screens(
                        beta0_calibrated, raster_info, cal_screens, cal_screens_raster_info, master_id
                    )
                    logging.info("...done.\n")

                elif self.proc_conf.processing_flags.DEM_flattening:
                    logging.info("FH: DEM flattening... ")
                    beta0_calibrated = apply_dem_flattening(
                        beta0_calibrated, kz, reference_height, master_id, raster_info
                    )
                    logging.info("...done.\n")

            except Exception as e:
                logging.error("FH: error during screen calibration or DEM flattening." + str(e), exc_info=True)
                raise

            ### compute mean incidence angle
            look_angle_rad = np.nanmean(off_nadir_angle_rad[master_id])  # 0.4886921905584123
            logging.info("FD: incidence angle used is {} [deg] \n".format(np.rad2deg(look_angle_rad)))

            ### ground notching
            try:
                logging.info("FD: ground contribute cancellation...:")

                beta0_notched = ground_cancellation(
                    beta0_calibrated,
                    kz,
                    self.proc_conf.ground_cancellation.multi_master_flag,
                    self.proc_conf.ground_cancellation.enhanced_forest_height,
                    self.proc_conf.ground_cancellation.equalization_flag,
                    raster_info.resolution_m_slant_rg,
                    off_nadir_angle_rad[master_id],
                    slope,
                )  # beta0_calibrated: 3 pol (nrg x Naz  x Nimmagini); beta0_notched: 3 pol (Nrg x Naz x 1 immagine)

                del beta0_calibrated, kz

            except Exception as e:
                logging.error("FD: error during ground cancellation: " + str(e), exc_info=True)
                raise
            logging.info("...done.\n")

            ### saving breakpoints
            if self.proc_conf.processing_flags.save_breakpoints:
                logging.info("FD: saving breakpoints (in slant range geometry) on " + breakpoints_output_folder)
                post_string = "_SR_" + unique_stack_id

                breakpoint_names = ["ground_cancelled" + post_string]

                save_breakpoints(breakpoints_output_folder, breakpoint_names, [beta0_notched])
                logging.info("...done.\n")

            # covariance estimation window size, it may be modified by an internal flag in case of air-plane geometry
            cov_est_window_size = self.proc_conf.change_detection_fd.product_resolution

            if self.proc_conf.processing_flags.multilook_heading_correction:
                _, heading_deg, _, _, _, _ = decode_unique_acquisition_id_string(unique_stack_id + "_BSL_00")

                cov_est_window_size = resolution_heading_correction(cov_est_window_size, heading_deg)

            logging.info("Covariance Estimation Window Size used: {} [m]".format(cov_est_window_size))

            try:
                logging.info("FD: disturbance covariance matrix computation...")
                (MPMB_covariance_sr, rg_vec_subs, az_vec_subs, subs_F_r, subs_F_a,) = main_covariance_estimation_SR(
                    deepcopy(beta0_notched),
                    cov_est_window_size,
                    raster_info.pixel_spacing_slant_rg,
                    raster_info.pixel_spacing_az,
                    look_angle_rad,
                    raster_info.carrier_frequency_hz,
                    raster_info.range_bandwidth_hz,
                )

                del beta0_notched

                MPMB_covariance_sr = covariance_matrix_mat2vec(
                    MPMB_covariance_sr
                )  # MPMB_covariance isreshaped from 3x3xNrgxNaz to 6xNrgxNaz

                # MPMB_covariance has 6 layers: saving each layer to numpy array and cancel the variable to free RAM space
                cov_layers_num, Nrg, Naz = MPMB_covariance_sr.shape  # cov_layers_num = 6
                MPMB_sr_temp_npy_fnames = []
                for layer_idx in np.arange(6):
                    MPMB_sr_temp_npy_fnames.append(
                        os.path.join(temp_output_folder_sr, "MPMB_covariance_layer_{}_6.npy".format(layer_idx + 1))
                    )
                    np.save(MPMB_sr_temp_npy_fnames[layer_idx], MPMB_covariance_sr[layer_idx, :, :])

                # del MPMB_covariance_sr
                logging.info("...done.\n")

            except Exception as e:
                logging.error("FD: error during covariance computation " + str(e), exc_info=True)
                raise

            ### Interpolate it over a regular lat lon grid (with grid data): generate a regular grid for the interpolation by using Max and Min lat lon from the ECEFGRID_LLH (make the vectors a bit longer), for the grid steps use the minimum steps from the ECEFGRID_LLH)
            logging.info(unique_stack_id + ": Geocoding data...")
            try:

                # initialize the geocoding
                min_spacing_m = min(
                    subs_F_a * raster_info.pixel_spacing_az, subs_F_r * raster_info.pixel_spacing_slant_rg
                )
                min_spacing_m = min(min_spacing_m, self.equi7_sampling)

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

                # geocode the MPMB_covariance
                logging.info(unique_stack_id + ": geocoding the covariance...")
                # geocode one layer at a time, loading from temporary numpy
                MPMB_gr_temp_npy_fnames = []
                for layer_idx in np.arange(cov_layers_num):
                    curr_geocoded_fname = os.path.join(
                        temp_output_folder_gr, "MPMB_covariance_layer_{}_6.npy".format(layer_idx + 1)
                    )

                    # load slant range MPMB from npy file geocode and save the geocoded on  a new npy file
                    MPMB_covariance_gr = geocoding(
                        np.load(MPMB_sr_temp_npy_fnames[layer_idx]),
                        lon_in,
                        lat_in,
                        lonMeshed_out,
                        latMeshed_out,
                        valid_values_mask,
                    )
                    np.save(curr_geocoded_fname, MPMB_covariance_gr)
                    MPMB_gr_temp_npy_fnames.extend([curr_geocoded_fname])

                del MPMB_covariance_gr

                logging.info("...done.\n")

                # geocode the incidence angle
                logging.info(unique_stack_id + ": geocoding the incidence angle...")
                # incidence angle = off nadir angle - slope
                inc_angle_ground_rad = geocoding(
                    (
                        off_nadir_angle_rad[master_id][rg_vec_subs, :][:, az_vec_subs]
                        - slope[rg_vec_subs, :][:, az_vec_subs]
                    ),
                    lon_in,
                    lat_in,
                    lonMeshed_out,
                    latMeshed_out,
                    valid_values_mask,
                )
                logging.info("...done.\n")

            except Exception as e:
                logging.error("FD: error during geocoding: " + str(e), exc_info=True)
                raise

            ### create GEOTIFF of all the layers (estimation, mask ):
            logging.info(unique_stack_id + ": formatting data to GEOTIFF...")

            try:

                cov_ground_fname = os.path.join(temp_output_folder_gr, "AverageCovariance.tif")
                inc_angle_ground_fname = os.path.join(temp_output_folder_gr, "IncidenceAngle.tif")

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

                # FD covariance matrix geotiff formatting: six layers in the same geotiff (multilayer)
                tiff_formatter(
                    MPMB_gr_temp_npy_fnames,
                    cov_ground_fname,
                    geotransform,
                    gdal_data_format=gdal.GDT_CFloat32,
                    multi_layers_tiff=True,
                    time_tag=str(time_tag_mjd_curr),
                )

                # incidence angle matrix geotiff formatting:
                tiff_formatter(
                    inc_angle_ground_rad,
                    inc_angle_ground_fname,
                    geotransform,
                    gdal_data_format=gdal.GDT_Float32,
                    multi_layers_tiff=False,
                    time_tag=str(time_tag_mjd_curr),
                )
            except Exception as e:
                logging.error("FD: error during GEOTIFF formatting: " + str(e), exc_info=True)
                raise

            del inc_angle_ground_rad

            ### formatting data to EQUI7 and writing the XML LookUp Tables for fast indexing
            logging.info(unique_stack_id + ": formatting into EQUI7 grid...")
            try:

                equi7_COV_parent_tempdir = os.path.join(temp_output_folder_e7, "AverageCovariance")

                # in general from here the Equi7 can output multiple tiles, which file names are stored in the output list ( wrapped here in a dict for the stack )
                equi7_cov_temp_fnames = image2equi7grid(
                    e7g,
                    cov_ground_fname,
                    equi7_COV_parent_tempdir,
                    gdal_path=gdal_path,
                    inband=None,
                    subgrid_ids=None,
                    accurate_boundary=False,
                    withtilenamesuffix=False,
                    resampling_type="bilinear",
                    tile_nodata=np.nan,
                )

                equi7_inc_parent_tempdir = os.path.join(temp_output_folder_e7, "IncidenceAngle")
                equi7_inc_temp_fnames = image2equi7grid(
                    e7g,
                    inc_angle_ground_fname,
                    equi7_inc_parent_tempdir,
                    gdal_path=gdal_path,
                    inband=None,
                    subgrid_ids=None,
                    accurate_boundary=False,
                    withtilenamesuffix=False,
                    resampling_type="bilinear",
                    tile_nodata=np.nan,
                )

            except Exception as e:
                logging.error("FD: error during EQUI7 formatting: " + str(e), exc_info=True)
                raise

            if not os.path.exists(equi7_cov_temp_fnames[0]):
                error_message = "EQUI7 grid has not been generated, output is absent "
                logging.error(error_message)
                raise
            else:
                data_driver = gdal.Open(equi7_cov_temp_fnames[0], GA_ReadOnly)
                projection_equi7 = data_driver.GetProjection()
                data_driver = None

            logging.info("...done.\n")

            if time_step_idx == 0:
                recomputed_fnf_folder = os.path.join(temp_output_folder, "input_fnf", "equi7")
                if os.path.exists(recomputed_fnf_folder):
                    # do not recompute fnf mask in output sampling if already present (because of precious cycles computation)
                    equi7_initial_fnf_mask_fnames = get_equi7_tiff_names(recomputed_fnf_folder)

                else:
                    # prepare the forest non forest mask
                    # it is loaded if equi7, or converted to equi7 if TANDEM-X
                    # it is a list containing all the loaded FNF-FTILES
                    fnf_format = check_fnf_folder_format(
                        self.proc_inputs.stack_based_processing.forest_mask_catalogue_folder
                    )
                    if fnf_format == "TANDEM-X":
                        logging.info("Initial Forest mask is in TANDEM-X format, converting to equi7...")

                        # conversion step 1: get the tiff names of current zone
                        equi7_initial_fnf_mask_fnames = fnf_tandemx_load_filter_equi7format(
                            self.proc_inputs.stack_based_processing.forest_mask_catalogue_folder,
                            e7g,
                            self.proc_conf.change_detection_fd.product_resolution,
                            temp_output_folder,
                            gdal_path,
                            geographic_boundaries,
                            self.time_tag_mjd_initial,
                        )

                    elif fnf_format == "EQUI7":
                        logging.info("Initial Forest mask reading and formatting...")

                        # re format the input mask in equi7 with the output resolution:
                        equi7_initial_fnf_mask_fnames = fnf_equi7_load_filter_equi7format(
                            self.proc_inputs.stack_based_processing.forest_mask_catalogue_folder,
                            e7g,
                            self.proc_conf.change_detection_fd.product_resolution,
                            temp_output_folder,
                            gdal_path,
                        )

            # check the equi7_fnf_mask_fnames
            check_equi7_mask_coverage(equi7_cov_temp_fnames, equi7_initial_fnf_mask_fnames)

            # equi7 zone is one contineltal zone
            equi7_zone_name = os.path.basename(os.path.dirname(os.path.dirname(equi7_cov_temp_fnames[0])))

            # current equi7 zone covariance product output folder
            equi7_zone_avg_cov_outdir = os.path.join(
                products_folder, "global_FD", "AverageCovariance", unique_stack_id, equi7_zone_name
            )
            equi7_zone_avg_cov_outdir_prev = os.path.join(
                products_folder, "global_FD", "AverageCovariance", unique_stack_id_prev, equi7_zone_name
            )

            # the disturbance computation is to be performed separately for each EQUI7 tile
            number_of_equi7_tiles = len(equi7_cov_temp_fnames)
            for equi7_tile_idx in np.arange(number_of_equi7_tiles):

                logging.info("FD: processing equi7 tile #{} of #{}".format(equi7_tile_idx + 1, number_of_equi7_tiles))

                # managing input/output paths and file names:
                equi7_tile_name = os.path.basename(
                    os.path.dirname(equi7_cov_temp_fnames[equi7_tile_idx])
                )  # current equi7 tile name inside the equi7 zone

                # load the mask for the equi7 tile if not already present (calling it "final" because it will be updated at each global cycle)
                if not equi7_tile_name in self.final_forest_mask_data.keys():

                    equi7_initial_fnf_mask_fname_curr = [
                        name for name in equi7_initial_fnf_mask_fnames if equi7_tile_name in name
                    ][0]

                    data_driver = gdal.Open(equi7_initial_fnf_mask_fname_curr, GA_ReadOnly)

                    self.final_forest_mask_data[equi7_tile_name] = data_driver.GetRasterBand(1).ReadAsArray()

                    time_tag_mjd_equi7_dict = data_driver.GetMetadata("TIFFTAG_DATETIME")
                    data_driver = None
                    if time_tag_mjd_equi7_dict:
                        time_tag_mjd_equi7 = PreciseDateTime().set_from_utc_string(time_tag_mjd_equi7_dict["time_tag"])

                        if time_tag_mjd_equi7 > self.time_tag_mjd_initial:
                            error_msg = 'Input FNF Mask (EQUI7) date of "{}" is bigger than input stack data minimum date of "{}"'.format(
                                str(time_tag_mjd_equi7), str(self.time_tag_mjd_initial)
                            )
                            logging.error(error_msg)
                            raise ValueError(error_msg)

                equi7_tile_avg_cov_outdir = os.path.join(
                    equi7_zone_avg_cov_outdir, equi7_tile_name
                )  # full path of equi7 tile output folder

                equi7_cov_tiff_name_curr = "AverageCovariance.tif"  # current equi7 tiff file name inside the equi7 tile

                equi7_avg_cov_out_tiff_name_prev = os.path.join(
                    equi7_zone_avg_cov_outdir_prev, equi7_tile_name, equi7_cov_tiff_name_curr
                )  # full path of equi7 tiff output file name at previous step (if present)

                if not os.path.exists(equi7_tile_avg_cov_outdir):
                    os.makedirs(equi7_tile_avg_cov_outdir)

                # compure number of looks (matrix Nr)
                data_driver = gdal.Open(equi7_inc_temp_fnames[equi7_tile_idx], GA_ReadOnly)
                inc_angle_equi7_rad = data_driver.GetRasterBand(1).ReadAsArray()
                mat_number_of_looks = np.floor(
                    (self.proc_conf.change_detection_fd.product_resolution ** 2)
                    / (raster_info.resolution_m_az * raster_info.resolution_m_slant_rg / np.sin(inc_angle_equi7_rad))
                )
                np.save(os.path.join(temp_output_folder, "mat_number_of_looks.npy"), mat_number_of_looks)

                # at first iteration we should load the current covariance matrix (Yn)
                # if it is not present, use the first MPMB_covariance Xi as Yn = Xi and pass to next stack
                if time_step_idx == 0:
                    prev_global_cycle_id = global_cycle_idx - 1
                    prev_global_cycle_id_str = str(prev_global_cycle_id)

                    equi7_avg_cov_out_tiff_name_prev = os.path.join(
                        self.proc_inputs.stack_based_processing.average_covariance_folder,
                        unique_stack_id_prev,
                        equi7_zone_name,
                        equi7_tile_name,
                        equi7_cov_tiff_name_curr,
                    )

                if time_step_idx == 0 and not os.path.exists(equi7_avg_cov_out_tiff_name_prev):

                    if prev_global_cycle_id < 0:
                        logging.info(
                            'FD first iteration, since first Global Cycle in input data is "0", cannot retrieve any previous Global Cycle Average Covariance'
                        )
                    else:
                        logging.info(
                            "FD first iteration, input Average Covariance tiff not present at: {}".format(
                                equi7_avg_cov_out_tiff_name_prev
                            )
                        )

                    logging.info(
                        "Saving current computed covariance as first Average Covariance (cannot compute disturbance at first iteration without input Average Covariance )"
                    )

                    # save the covariance to output and continue to next global cycle iteration
                    shutil.move(equi7_cov_temp_fnames[equi7_tile_idx], equi7_tile_avg_cov_outdir)

                    continue

                elif time_step_idx == 0 and os.path.exists(equi7_avg_cov_out_tiff_name_prev):

                    logging.info(
                        "FD: first iteration, loading Average Covariance from of previous Global Cycle {}".format(
                            prev_global_cycle_id_str
                        )
                    )

                    # the starting Yn is an input:
                    data_driver = gdal.Open(equi7_avg_cov_out_tiff_name_prev, GA_ReadOnly)
                    for layer_idx in np.arange(6):
                        # saved matrix is 6xNrgxNaz, instead Y_N_matrix is 3x3xNrgxNaz matrix
                        if layer_idx == 0:
                            layer_read = data_driver.GetRasterBand(int(layer_idx + 1)).ReadAsArray()
                            Nrg_equi7, Naz_equi7 = layer_read.shape
                            Y_N_6x_vec = np.zeros((6, Nrg_equi7, Naz_equi7), dtype=np.complex64)
                            Y_N_6x_vec[layer_idx, :, :] = layer_read
                        else:
                            Y_N_6x_vec[layer_idx, :, :] = data_driver.GetRasterBand(int(layer_idx + 1)).ReadAsArray()

                    data_driver = None

                logging.info('FD, global cycle "GC_' + global_cycle_idx_str + '" : computing distubance...')

                if time_step_idx != 0:
                    # loading cumulate covariance at previous step:
                    data_driver = gdal.Open(equi7_avg_cov_out_tiff_name_prev, GA_ReadOnly)

                    for layer_idx in np.arange(6):

                        if layer_idx == 0:
                            layer_read = data_driver.GetRasterBand(int(layer_idx + 1)).ReadAsArray()
                            Nrg_equi7, Naz_equi7 = layer_read.shape
                            Y_N_6x_vec = np.zeros((6, Nrg_equi7, Naz_equi7), dtype=np.complex64)
                            Y_N_6x_vec[layer_idx, :, :] = layer_read
                        else:
                            Y_N_6x_vec[layer_idx, :, :] = data_driver.GetRasterBand(int(layer_idx + 1)).ReadAsArray()

                    del layer_read
                    data_driver = None

                # loading covariance at current step:
                data_driver = gdal.Open(equi7_cov_temp_fnames[equi7_tile_idx], GA_ReadOnly)

                for layer_idx in np.arange(6):

                    if layer_idx == 0:
                        layer_read = data_driver.GetRasterBand(int(layer_idx + 1)).ReadAsArray()
                        X_i_6x_vec = np.zeros((6, Nrg_equi7, Naz_equi7), dtype=np.complex64)
                    else:
                        X_i_6x_vec[layer_idx, :, :] = data_driver.GetRasterBand(int(layer_idx + 1)).ReadAsArray()

                del layer_read
                data_driver = None

                j_idx = 1
                disturbance_matrix = np.zeros((Nrg_equi7, Naz_equi7), dtype=np.bool)
                disturbance_prob_matrix = np.zeros((Nrg_equi7, Naz_equi7), dtype=np.bool)

                Nrg_equi7_str = str(Nrg_equi7)
                for rg_idx in np.arange(Nrg_equi7):

                    if not np.remainder(rg_idx, 50):
                        logging.info("    change detection step {} of ".format(rg_idx) + Nrg_equi7_str)

                    for az_idx in np.arange(Naz_equi7):

                        if np.any(np.isnan(X_i_6x_vec[:, rg_idx, az_idx])):
                            continue

                        X_i = covariance_matrix_vec2mat(X_i_6x_vec[:, rg_idx, az_idx])  # X_i is 3x3

                        Y_N = covariance_matrix_vec2mat(Y_N_6x_vec[:, rg_idx, az_idx])  # Y_N is 3x3
                        cov_number_of_looks = mat_number_of_looks[rg_idx, az_idx]

                        (
                            Y_N,
                            j_idx,
                            disturbance_matrix[rg_idx, az_idx],
                            disturbance_prob_matrix[rg_idx, az_idx],
                            _,
                        ) = alg_wishart_SU(
                            X_i,
                            Y_N,
                            j_idx,
                            number_of_pols,
                            cov_number_of_looks,
                            self.proc_conf.change_detection_fd.confidence_level,
                            False,
                        )

                        # update the matrix for the next stack
                        Y_N_6x_vec[:, rg_idx, az_idx] = covariance_matrix_mat2vec(Y_N)

                logging.info("...done.\n")

                logging.info("updating fnf mask for current global cycle")
                self.final_forest_mask_data[equi7_tile_name][disturbance_matrix] = False

                logging.info("saving updated average covariance for current global cycle")

                # SAVE Equi7 Average Covariance
                equi7_initial_fnf_mask_fname_curr = [
                    name for name in equi7_initial_fnf_mask_fnames if equi7_tile_name in name
                ][0]
                data_driver = gdal.Open(equi7_initial_fnf_mask_fname_curr, GA_ReadOnly)
                geotransform_out = data_driver.GetGeoTransform()
                data_driver = None

                equi7_COV_out_fname = os.path.join(equi7_tile_avg_cov_outdir, equi7_cov_tiff_name_curr)
                tiff_formatter(
                    [
                        Y_N_6x_vec[0, :, :],
                        Y_N_6x_vec[1, :, :],
                        Y_N_6x_vec[2, :, :],
                        Y_N_6x_vec[3, :, :],
                        Y_N_6x_vec[4, :, :],
                        Y_N_6x_vec[5, :, :],
                    ],
                    equi7_COV_out_fname,
                    geotransform_out,
                    gdal_data_format=gdal.GDT_Float32,
                    projection=projection_equi7,
                    multi_layers_tiff=True,
                    time_tag=str(time_tag_mjd_curr),
                )

                logging.info(
                    "FD: equi7 tile {} of {}".format(equi7_tile_idx, len(equi7_cov_temp_fnames))
                    + " average covariange processing done.\n"
                )

            # keep track of the unique_stack_id for next iteration
            unique_stack_id_prev = unique_stack_id

            logging.info(
                'FD: disturbance computation for global cycle "GC_'
                + global_cycle_idx_str
                + '" of geometry '
                + self.nominal_geometry_stack_id
                + " done.\n"
            )

        ################## CYCLES LOOP END ##############################

        # save the final mask (one for each equi7 tile) and compute the dsturbance (one for each equi7 tile too)
        time_tag_mjd_final = time_tag_mjd_curr
        equi7_zone_disturbance_outdir = os.path.join(
            products_folder, "global_FD", "Disturbance", self.nominal_geometry_stack_id, equi7_zone_name
        )
        equi7_zone_fnf_mask_outdir = os.path.join(
            products_folder, "global_FD", "FNF_Mask", self.nominal_geometry_stack_id, equi7_zone_name
        )

        for equi7_tile_idx in np.arange(number_of_equi7_tiles):

            equi7_tile_name = os.path.basename(os.path.dirname(equi7_cov_temp_fnames[equi7_tile_idx]))

            equi7_tile_disturbance_outdir = os.path.join(equi7_zone_disturbance_outdir, equi7_tile_name)

            equi7_fnf_mask_tile_outdir = os.path.join(equi7_zone_fnf_mask_outdir, equi7_tile_name)
            equi7_disturbance_tiff_name = "disturbance.tif"
            equi7_fnf_mask_tiff_name = "fnf.tif"

            equi7_disturbance_out_fname = os.path.join(equi7_tile_disturbance_outdir, equi7_disturbance_tiff_name)
            equi7_fnf_mask_outfname = os.path.join(equi7_fnf_mask_tile_outdir, equi7_fnf_mask_tiff_name)

            if not os.path.exists(equi7_tile_disturbance_outdir):
                os.makedirs(equi7_tile_disturbance_outdir)
            if not os.path.exists(equi7_fnf_mask_tile_outdir):
                os.makedirs(equi7_fnf_mask_tile_outdir)

            logging.info(
                "FD: final step, computation of disturbance product, for equi7 tile #{} of #{}".format(
                    equi7_tile_idx + 1, number_of_equi7_tiles
                )
            )

            equi7_initial_fnf_mask_fname_curr = [
                name for name in equi7_initial_fnf_mask_fnames if equi7_tile_name in name
            ][0]
            data_driver = gdal.Open(equi7_initial_fnf_mask_fname_curr, GA_ReadOnly)
            initial_forest_mask_curr = data_driver.GetRasterBand(1).ReadAsArray()
            geotransform_out = data_driver.GetGeoTransform()
            data_driver = None

            final_disturbance_matrix_curr = np.logical_xor(
                initial_forest_mask_curr, self.final_forest_mask_data[equi7_tile_name]
            )

            tiff_formatter(
                final_disturbance_matrix_curr,
                equi7_disturbance_out_fname,
                geotransform_out,
                gdal_data_format=gdal.GDT_Float32,
                projection=projection_equi7,
                time_tag=str(time_tag_mjd_final),
            )

            logging.info(
                "FD: final step,saving the final fnf mask for equi7 tile #{} of #{}".format(
                    equi7_tile_idx, number_of_equi7_tiles
                )
            )
            tiff_formatter(
                self.final_forest_mask_data[equi7_tile_name],
                equi7_fnf_mask_outfname,
                geotransform_out,
                gdal_data_format=gdal.GDT_Float32,
                projection=projection_equi7,
                time_tag=str(time_tag_mjd_final),
            )

        # return the final_forest_mask_data to be updated in the next stack-cycle
        return self.final_forest_mask_data, temp_output_folder
