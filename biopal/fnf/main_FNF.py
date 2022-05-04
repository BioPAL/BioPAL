# SPDX-FileCopyrightText: BioPAL <biopal@esa.int>
# SPDX-License-Identifier: MIT

import os
import numpy as np
import logging
from scipy.signal import convolve2d
import shutil
from osgeo import gdal
from osgeo.gdalconst import GA_ReadOnly
from equi7grid.equi7grid import Equi7Grid
from equi7grid.image2equi7grid import image2equi7grid


from biopal.data_operations.data_operations import (
    read_and_oversample_data,
    read_and_oversample_aux_data,
    apply_dem_flattening,
    mosaiking,
)
from biopal.utility.utility_functions import (
    Task,
    set_gdal_paths,
    choose_equi7_sampling,
    check_if_path_exists,
    check_if_geometry_auxiliaries_are_present,
    resolution_heading_correction,
    decode_unique_acquisition_id_string,
    evaluate_estimation_quality_matrix,
    save_breakpoints,
    collect_stacks_to_be_merged,
)
from biopal.geocoding.geocoding import (
    geocoding,
    geocoding_init,
)
from biopal.io.xml_io import (
    parse_input_file,
    parse_configuration_file,
    write_input_file,
    core_processing_fnf,
)
from biopal.statistics.utility_statistics import (
#     build_filtering_matrix,
#     Covariance2D2Correlation2D,
#    main_correlation_estimation_SR,
    main_covariance_estimation_SR,
#     covariance_matrix_vec2mat,
)
from biopal.io.data_io import tiff_formatter
from biopal.screen_calibration.screen_calibration import apply_calibration_screens
from biopal.geometry.utility_geometry import compute_and_oversample_geometry_auxiliaries
from biopal.fnf.processing_FNF import estimate_pixel_coefficients
from biopal.fnf.processing_FNF import get_centered_logreg_measures
from biopal.fnf.processing_FNF import logistic_regression




def pforest_masking_and_merging(data_equi7_fnames, mask_equi7_fnames, stacks_to_merge_dict):

    # for each stack:
    # 1) check if it is alone, of if there is a couple ASC+DES stacks
    # 2) if alone, just apply the mask
    # 3) if it is not alone, merge it togheter with its companion data, also applying the mask

    # the asc_des_string is formatted as it follow:
    # XXX:SELF,YYY:STACK_ID
    # where:
    #   XXX can be 'ASC' or 'DES'
    #   YYY can be 'ASC' or 'DES' (if XXX='ASC', than YYY='DES' and vicecersa )
    #   'SELF' is a constant string: the asc_des_string is contained in a dictionary, 'SELF' is the dictionary key
    #   STACK_ID can be 'N.A' or can be the name of a stack_id

    # each input in data_equi7_fnames is a tiff with two layers
    # first later is data
    data_layer_index = 1
    # second layer is quality
    quality_layer_index = 2
    # third layer is FNF mask (0: NoForest, 1: Forest, nan: No valid value)
    fnf_layer_index = 3

    temp_path = next(iter(data_equi7_fnames.values()))[0]
    idx = temp_path.find("equi7")
    merging_folder = os.path.join(temp_path[: idx + len("equi7")], "merged")

    merged_data_fnames = []
    for (current_merging_id, unique_stack_ids_to_merge_list) in stacks_to_merge_dict.items():

        # prepare the empty dictionary for equi7 tile names
        equi7_tile_names = []

        # retrive the names of all the used "EQUI7 tiles"
        for curr_equi7_name in data_equi7_fnames[unique_stack_ids_to_merge_list[0]]:
            curr_equi7_tile_name = os.path.basename(os.path.dirname(curr_equi7_name))
            if not curr_equi7_tile_name in equi7_tile_names:
                equi7_tile_names.append(curr_equi7_tile_name)

        first_stack_data_fname = data_equi7_fnames[unique_stack_ids_to_merge_list[0]][0]
        data_tiff_name_in = os.path.basename(first_stack_data_fname)

        data_tiff_name_out = data_tiff_name_in.replace(unique_stack_ids_to_merge_list[0],
                                                       current_merging_id)

        p1 = os.path.dirname(os.path.dirname(first_stack_data_fname))
        equi7_subgridname = os.path.basename(p1)

        for equi7_tile_name in equi7_tile_names:

            out_data_fname = os.path.join(
                merging_folder, current_merging_id, equi7_subgridname,
                equi7_tile_name, data_tiff_name_out,
            )

            os.makedirs(os.path.dirname(out_data_fname))

            for idx, unique_stack_id in enumerate(unique_stack_ids_to_merge_list):
                # input and output paths

                curr_data_fname = [name for name in data_equi7_fnames[unique_stack_id] if equi7_tile_name in name][0]
                curr_mask_fname = [name for name in mask_equi7_fnames[unique_stack_id] if equi7_tile_name in name][0]

                # Load DATA
                data_driver = gdal.Open(curr_data_fname, GA_ReadOnly)
                data_curr = data_driver.GetRasterBand(data_layer_index).ReadAsArray()
                projection = data_driver.GetProjection()
                geotransform = data_driver.GetGeoTransform()
                # Load quality layer
                quality_data_curr = data_driver.GetRasterBand(quality_layer_index).ReadAsArray()
                data_driver = None

                # Load MASK
                mask_driver = gdal.Open(curr_mask_fname, GA_ReadOnly)
                mask_curr = (mask_driver.GetRasterBand(1).ReadAsArray()).astype("bool")
                mask_driver = None

                if idx == 0:
                    Pforest_sum_data = np.zeros_like(data_curr)
                    quality_sum_data = np.zeros_like(quality_data_curr)
                    N_values = np.zeros_like(Pforest_sum_data)
                    Pforest_sum_data[mask_curr] = data_curr[mask_curr]
                    quality_sum_data[mask_curr] = quality_data_curr[mask_curr]
                    N_values[mask_curr] += 1

                else:

                    # data_prev = Pforest_sum_data
                    # quality_data_prev = quality_sum_data
                    # del Pforest_out_data, quality_out_data

                    # 1) Add the current value to previous:
                    Pforest_sum_data[mask_curr] += data_curr[mask_curr]
                    quality_sum_data[mask_curr] += quality_data_curr[mask_curr]
                    # 2) increase the number of valid points for each pixel
                    N_values[mask_curr] += 1

            # Make the mean and generate global mask
            # Total mask with the pixels having at least 1 value
            total_mask = (N_values > 0)
            Pforest_avg = np.zeros_like(Pforest_sum_data)
            quality_avg = np.zeros_like(quality_sum_data)
            # Perform the mean
            Pforest_avg[total_mask] = Pforest_sum_data[total_mask] / N_values[total_mask]
            quality_avg[total_mask] = quality_sum_data[total_mask] / N_values[total_mask]
            # Set values without data (N_values == 0) to np.nan
            Pforest_avg[~total_mask] = np.nan
            quality_avg[~total_mask] = np.nan
            # Generate FNF mask by thresholding Pforest
            FNF_mask = (Pforest_avg > 0.5).astype(np.float32)
            FNF_mask[~total_mask] = np.nan
            
            # save to file:
            Nx, Ny = Pforest_avg.shape
            driver = gdal.GetDriverByName("GTiff")
            outdata = driver.Create(out_data_fname, Ny, Nx, 3, gdal.GDT_Float32)
            outdata.SetGeoTransform(geotransform)  ##sets same geotransform as input
            outdata.SetProjection(projection)  ##sets same projection as input
            outdata.GetRasterBand(data_layer_index).WriteArray(Pforest_avg)
            outdata.GetRasterBand(quality_layer_index).WriteArray(quality_avg)
            outdata.GetRasterBand(fnf_layer_index).WriteArray(FNF_mask)
            outdata.FlushCache()  ##saves to disk!!
            outdata = None

        merged_data_fnames.append(out_data_fname)
        logging.info("    ...done.")

    return merged_data_fnames, merging_folder



class ForestNonForestMap(Task):
    """FNF main APP ForestNonForestMap

    run this APP to execute the complete Forest Non-Forest map processing chain.
    
    ForestNonForestMap is composed by two sub APPS automatically called in sequence  when standard launch is performed:
    StackBasedProcessingFNF -> CoreProcessingFNF
   
    Refer to dataset_query, StackBasedProcessingFNF and CoreProcessingFNF documentation for step by step run.

    Attributes
    ----------
    configuration_file : str
        path of the Configuration_File.xml file

    Methods
    -------
    run( input_file_path )
        run the ForestNonForestMap processing
    name : str
        name of the APP
        
    See Also
    --------
    biopal.dataset_query.dataset_query.dataset_query: it's the APP to be called before this APP
    StackBasedProcessingFNF : it's the first of the two sub-APPs called by ForestNonForestMap preocessor
    CoreProcessingFNF : it's  the second of the two sub-APPs called by ForestNonForestMap processor

    Examples
    --------
    Manual FNF chain execution
    
    >>> from biopal.dataset_query.dataset_query import dataset_query
    >>> from biopal.fnf.main_FNF import ForestNonForestMap
    >>> dq_obj = dataset_query()
    >>> input_file_up = dq_obj.run( input_file )
    >>> chain_obj = ForestNonForestMap( configuration_file )
    >>> chain_obj.run( input_file_up )

    - input_file: path of the BioPAL input file
    - input_file_up: same of input_file with also the "stack_based_processing" section
    - configuration_file: path of the BioPAL configuration file
    """

    def __init__(self, configuration_file):

        super().__init__(configuration_file)

    def _run(self, input_file):

        # Main APP #1: Stack Based Processing
        stack_based_processing_obj = StackBasedProcessingFNF(self.configuration_file)

        # Run Main APP #1: Stack Based Processing
        input_file_updated = stack_based_processing_obj.run(input_file)

        # Main APP #2: Core Processing
        fnf_processing_obj = CoreProcessingFNF(self.configuration_file)

        # Run Main APP #2: FNF Core Processing
        fnf_processing_obj.run(input_file_updated)


class StackBasedProcessingFNF(Task):
    """StackBasedProcessingFNF APP 
    
    StackBasedProcessingFNF APP is the first of the two sub-APPs called by ForestNonForestMap processor.
    
    It performs the stack-based probability of forest estimation with logistic
    regression.

    Attributes
    ----------
    configuration_file : str
        path of the Configuration_File.xml file
        
    Methods
    -------
    run( input_file_path )
        run the StackBasedProcessingFNF APP
    name : str
        name of the APP
        
    See Also
    --------
    biopal.dataset_query.dataset_query.dataset_query : it's the first APP to be called in the manual sequence
    CoreProcessingFNF : it's the core APP that follows this APP in the call sequence

    Examples
    --------
    Manual FNF chain execution
        
    >>> from biopal.dataset_query.dataset_query import dataset_query
    >>> from biopal.fnf.main_FNF import StackBasedProcessingFNF
    >>> dq_obj = dataset_query()
    >>> input_file_up1 = dq_obj.run( input_file )
    >>> sbp_obj = StackBasedProcessingFNF( config_file )
    >>> input_file_up2 = sbp_obj.run( input_file_up1 )
    >>> fnfcp_obj = CoreProcessingFNF( config_file )
    >>> fnfcp_obj.run( input_file_up2 )

    - input_file: path of the BioPAL input file
    - input_file_up1: same of input_file with also the "stack_based_processing" section   
    - input_file_up2: same of input_file_up1 with also the "core_processing_fnf" section
    - config_file: path of the BioPAL configuration file
    """

    def __init__(self, configuration_file):
        super().__init__(configuration_file)

    def _run(self, input_file):

        ########################## INITIAL STEPS ##############################

        logging.info("FNF: Reading chains configuration file")
        check_if_path_exists(self.configuration_file, "FILE")
        conf_params_obj = parse_configuration_file(self.configuration_file)
        input_params_obj = parse_input_file(input_file)

        ### managing output folders:
        products_folder = os.path.join(input_params_obj.output_specification.output_folder, "Products")
        if conf_params_obj.processing_flags.save_breakpoints:
            breakpoints_output_folder = os.path.join(products_folder, "breakpoints")
            logging.info("FNF: Breakpoints will be saved into: " + breakpoints_output_folder)
            os.makedirs(breakpoints_output_folder)

        temp_output_folder = os.path.join(products_folder, "temp")
        logging.info("FNF: Temporary data folder:" + temp_output_folder + "\n")
        os.makedirs(temp_output_folder)

        equi7_sampling = choose_equi7_sampling(
            conf_params_obj.estimate_fnf.product_resolution,
            input_params_obj.output_specification.geographic_grid_sampling,
        )
        e7g = Equi7Grid(equi7_sampling)
        logging.info("    EQUI7 Grid sampling used: {}".format(equi7_sampling))

        # get needed parameters from input and configuration files
        gdal_path, _ = set_gdal_paths(conf_params_obj.gdal.gdal_path, conf_params_obj.gdal.gdal_environment_path)

        ########################## INITIAL STEPS END #############################

        pforest_equi7_fnames = {}
        mask_equi7_fnames = {}
        ########################## STACK BASED STEPS ##############################
        for unique_stack_id, acquisitions_pf_names in input_params_obj.stack_based_processing.stack_composition.items():

            # make temporary sub-directories
            temp_output_folder_gr = os.path.join(temp_output_folder, "geocoded", unique_stack_id)
            temp_output_folder_e7 = os.path.join(temp_output_folder, "equi7", unique_stack_id)
            os.makedirs(temp_output_folder_gr)
            os.makedirs(temp_output_folder_e7)

            ### load data ( and oversample if requested and if needed )
            try:
                logging.info("FNF: Data loading for stack " + unique_stack_id + "; this may take a while:")

                (beta0_calibrated, master_id, raster_info, raster_info_orig,) = read_and_oversample_data(
                    input_params_obj.dataset_query.L1C_repository,
                    acquisitions_pf_names,
                    conf_params_obj.processing_flags.enable_resampling,
                )
            except Exception as e:
                logging.error("FNF: error during input data reading: " + str(e), exc_info=True)
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
                    read_dist=False,
                )

                if conf_params_obj.processing_flags.compute_geometry or not geometry_aux_are_present:

                    # messages for the log:
                    if conf_params_obj.processing_flags.compute_geometry:
                        logging.info("FNF: calling geometry library for stack " + unique_stack_id + "\n")
                        if geometry_aux_are_present:
                            logging.warning("    geometry auxiliaries will be overwritten for stack " + unique_stack_id)
                        else:
                            logging.info("\n")
                    else:
                        logging.warning(
                            'FNF: calling geometry library since AuxiliaryProductsFolder "Geometry" is empty or not complete \n'
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
                        comp_ref_h=read_ref_h,
                        comp_dist=False,
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
                        comp_dist=False,
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

                    logging.info("FNF: ...geometry auxiliaries computation done.")

                else:

                    logging.info(
                        "FNF: geometry auxiliaries are provided from user, so they are now being loaded and not computed, for stack "
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
                        read_dist=False,
                    )
                    logging.info("FNF: ...geometry auxiliaries loading done.")

                # read the rest of auxiliaries which are notpart of the geometry library:
                if read_cal_screens:

                    logging.warning("FNF: loading calibration screens \n")

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
                    "FNF: error during auxiliary data computation and/or loading: " + str(e), exc_info=True,
                )
                raise

            ### Screen calibration (ground steering)
            try:
                if conf_params_obj.processing_flags.apply_calibration_screen:
                    logging.info("FNF: applying calibration screen...")
                    beta0_calibrated = apply_calibration_screens(
                        beta0_calibrated, raster_info, cal_screens, cal_screens_raster_info, master_id,
                    )
                    logging.info("...done.\n")

                elif conf_params_obj.processing_flags.DEM_flattening:
                    logging.info("FNF: DEM flattening... ")
                    beta0_calibrated = apply_dem_flattening(
                        beta0_calibrated, kz, reference_height, master_id, raster_info
                    )
                    logging.info("...done.\n")

            except Exception as e:
                logging.error(
                    "FNF: error during screen calibration or DEM flattening." + str(e), exc_info=True,
                )
                raise

            ### compute mean incidence angle
            look_angle_rad = np.nanmean(off_nadir_angle_rad[master_id])  # 0.4886921905584123
            logging.info("FNF: incidence angle used is {} [deg] \n".format(np.rad2deg(look_angle_rad)))

            ### mean of off nadir and slope over final resolution
            windtm_x = np.int(
                np.round(conf_params_obj.estimate_fnf.product_resolution / raster_info.pixel_spacing_az / 2) * 2 + 1
            )
            windtm_y = np.int(
                np.round(
                    conf_params_obj.estimate_fnf.product_resolution
                    / (raster_info.pixel_spacing_slant_rg / np.sin(look_angle_rad))
                    / 2
                )
                * 2
                + 1
            )

            logging.info("FNF: computing mean of off nadir based on final product resolution...")

            # only mster is used in the code, do not convolve all the stacs, it is useless
            off_nadir_angle_rad[master_id] = convolve2d(
                off_nadir_angle_rad[master_id], np.ones((windtm_y, windtm_x)) / windtm_y / windtm_x, mode="same",
            )
            logging.info("...done.")

            logging.info("FNF: computing mean of slope based on final product resolution...")

            slope = convolve2d(slope, np.ones((windtm_y, windtm_x)) / windtm_y / windtm_x, mode="same")
            logging.info("...done.")

            # covariance estimation window size, it may be modified by an internal flag in case of air-plane geometry
            cov_est_window_size = conf_params_obj.estimate_fnf.product_resolution

            if conf_params_obj.processing_flags.multilook_heading_correction:
                _, heading_deg, _, _, _, _ = decode_unique_acquisition_id_string(unique_stack_id + "_BSL_00")

                cov_est_window_size = resolution_heading_correction(cov_est_window_size, heading_deg)

            logging.info("Covariance Estimation Window Size used: {} [m]".format(cov_est_window_size))

            ### Pforest estimation
            logging.info("FNF: " + unique_stack_id + ": performing probability of forest estimation...")
            try:

                (prob_forest, rg_vec_subs, az_vec_subs, subs_F_r, subs_F_a,) = EstimateProbForest(
                    beta0_calibrated,
                    cov_est_window_size,
                    raster_info.pixel_spacing_slant_rg,
                    raster_info.pixel_spacing_az,
                    look_angle_rad,
                    raster_info.carrier_frequency_hz,
                    raster_info.range_bandwidth_hz,
                    kz,
                    conf_params_obj.estimate_fnf.logreg_coeffs_fname,
                )

            except Exception as e:
                logging.error("FNF: error during probability of forest estimation: " + str(e), exc_info=True)
                raise
            del beta0_calibrated

            ### Placemark for the quality estimation to be defined
            try:
                quality_layer_sr = evaluate_estimation_quality_matrix(prob_forest.shape)

            except Exception as e:
                logging.error("FNF: error during estimation quality layer evaluation: " + str(e), exc_info=True)
                raise

            ### Interpolate it over a regular lat lon grid (with grid data): generate a regular grid for the interpolation by using Max and Min lat lon from the ECEFGRID_LLH (make the vectors a bit longer), for the grid steps use the minimum steps from the ECEFGRID_LLH)
            logging.info(unique_stack_id + ": Geocoding data...")
            try:

                # initialize the geocoding
                min_spacing_m = min(
                    subs_F_a * raster_info.pixel_spacing_az, subs_F_r * raster_info.pixel_spacing_slant_rg,
                )
                min_spacing_m = min(min_spacing_m, equi7_sampling)

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

                # geocode the Probability of Forest
                logging.info(unique_stack_id + ": geocoding the probability of forest...")
                pforest_ground = geocoding(
                    prob_forest, lon_in, lat_in, lonMeshed_out, latMeshed_out, valid_values_mask,
                )
                logging.info("...done.\n")

                # geocode the Probability of Forest quality layer
                logging.info(unique_stack_id + ": geocoding the estimated Pforest quality layer...")
                quality_layer_ground = geocoding(
                    quality_layer_sr, lon_in, lat_in, lonMeshed_out, latMeshed_out, valid_values_mask,
                )
                logging.info("...done.\n")

            except Exception as e:
                logging.error("FNF: error during geocoding: " + str(e), exc_info=True)
                raise

            ### saving breakpoints
            if conf_params_obj.processing_flags.save_breakpoints:
                logging.info("FNF: saving main results (in slant range geometry) on " + breakpoints_output_folder)
                post_string = "_SR_" + unique_stack_id

                breakpoint_names = ["probability_forest" + post_string]

                save_breakpoints(breakpoints_output_folder, breakpoint_names, [prob_forest])
                logging.info("...done.\n")
            del prob_forest

            ### creating mask to exclude estimation not valid values:
            # Currently it is always True
            condition_curr = pforest_ground >= 0.0
            
            estimation_mask_ground = np.where(condition_curr, True, False)
            # Mask points without data
            estimation_mask_ground[np.isnan(pforest_ground)] = False
            # also casted to float, for incapsulation in a geotiff
            estimation_mask_ground.astype(float)

            ### create GEOTIFF of all the layers (estimation, mask ):
            logging.info(unique_stack_id + ": formatting data to GEOTIFF...")
            try:

                pforest_ground_fname = os.path.join(temp_output_folder_gr, "Pforest.tif")
                valid_values_mask_ground_fname = os.path.join(temp_output_folder_gr, "Valid_mask.tif")

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

                # probability of forest geotiff formatting
                tiff_formatter(
                    [pforest_ground, quality_layer_ground],
                    pforest_ground_fname,
                    geotransform,
                    gdal_data_format=gdal.GDT_Float32,
                    multi_layers_tiff=True,
                )

                # FNF temporary mask geotiff formatting
                tiff_formatter(
                    estimation_mask_ground,
                    valid_values_mask_ground_fname,
                    geotransform,
                    gdal_data_format=gdal.GDT_Float32,
                )

                logging.info("...done.\n")

            except Exception as e:
                logging.error("FNF: error during GEOTIFF formatting: " + str(e), exc_info=True)
                raise

            ### formatting data to EQUI7
            logging.info(unique_stack_id + ": formatting into EQUI7 grid...")
            try:

                equi7_pforest_outdir = os.path.join(temp_output_folder_e7, "p_forest")
                equi7_mask_outdir = os.path.join(temp_output_folder_e7, "p_forest_mask")

                # in general from here the Equi7 can output multiple tiles, which file names are stored in the output list ( wrapped here in a dict for the stack )
                pforest_equi7_fnames[unique_stack_id] = image2equi7grid(
                    e7g,
                    pforest_ground_fname,
                    equi7_pforest_outdir,
                    gdal_path=gdal_path,
                    inband=None,
                    subgrid_ids=None,
                    accurate_boundary=False,
                    withtilenamesuffix=False,
                    resampling_type="bilinear",
                    tile_nodata=np.nan,
                )
                mask_equi7_fnames[unique_stack_id] = image2equi7grid(
                    e7g,
                    valid_values_mask_ground_fname,
                    equi7_mask_outdir,
                    gdal_path=gdal_path,
                    inband=None,
                    subgrid_ids=None,
                    accurate_boundary=False,
                    withtilenamesuffix=False,
                    resampling_type="bilinear",
                    tile_nodata=np.float(0),
                )

            except Exception as e:
                logging.error("FNF: error during EQUI7 formatting: " + str(e), exc_info=True)
                raise

            if not os.path.exists(pforest_equi7_fnames[unique_stack_id][0]):
                error_message = "EQUI7 grid has not been generated, output is absent "
                logging.error(error_message)
                raise

            logging.info("...done.\n")

        # write the input file with the sections needed by the Core Processing FNF APP:
        out_input_file_xml = os.path.join(
            os.path.dirname(input_params_obj.output_specification.output_folder), "Input_File_CoreProcessingFNF.xml"
        )
        input_params_obj.core_processing_fnf = fill_core_processing_fnf_obj(
            input_params_obj, pforest_equi7_fnames, mask_equi7_fnames,
        )
        write_input_file(input_params_obj, out_input_file_xml)

        ######################## STACK BASED STEPS END. ###########################
        return out_input_file_xml

def EstimateProbForest(
    data_stack,
    cov_est_window_size,
    pixel_spacing_slant_rg,
    pixel_spacing_az,
    incidence_angle_rad,
    carrier_frequency_hz,
    range_bandwidth_hz,
    kz_stack,
    logreg_coeffs_fname,
):


    # data_stack is a dictionary of two nested dictionaries composed as:
    # data_stack[ acquisition_name ][ polarization ]

    num_acq = len(data_stack)
    # acq_names = list(data_stack.keys())
    # first_acq_dict = data_stack[acq_names[0]]
    # pol_names = list(first_acq_dict.keys())
    # num_pols = len(pol_names)
    # Nrg, Naz = first_acq_dict[pol_names[0]].shape

    # Covariance estimation
    (MPMB_correlation, rg_vec_subs, az_vec_subs, subs_F_r, subs_F_a,) = main_covariance_estimation_SR(
        data_stack,
        cov_est_window_size,
        pixel_spacing_slant_rg,
        pixel_spacing_az,
        incidence_angle_rad,
        carrier_frequency_hz,
        range_bandwidth_hz,
    )

    # Nrg_subs = rg_vec_subs.size
    # Naz_subs = az_vec_subs.size
    
    # Covariance matrix indices on last axes, as expected by FNF core functions
    MPMB_correlation = np.transpose(MPMB_correlation, axes=(2,3,0,1))
    
    # Build kz array
    kz = np.zeros(MPMB_correlation.shape[:2] + (num_acq,) )
    for b_idx, stack_curr in enumerate(kz_stack.values()):
        kz[..., b_idx] = stack_curr[
            tuple(np.meshgrid(rg_vec_subs, az_vec_subs, indexing='ij'))]
    
    # Load coefficients table
    logging.debug("Loading logistic regresion coefficientes file '%s'..." %
                 logreg_coeffs_fname)
    with np.load(logreg_coeffs_fname) as f:
        log_reg_maxkz = f['log_reg_maxkz']
        log_reg_coeffs = f['log_reg_coeffs']
    

    logging.debug("Obtaining logistic regression features...")
    sv1, coh_std, pi, p12coh, hhvvcoh, trimg = get_centered_logreg_measures(MPMB_correlation)
    
    logging.debug("Computing coefficients for each pixel...")
    coefs = estimate_pixel_coefficients(kz, log_reg_coeffs, log_reg_maxkz)
    
    logging.debug("Generating feature vector...")
    x_full = np.stack([sv1, coh_std, pi, p12coh, hhvvcoh, trimg], axis=-1)
    logging.debug("... and calling logistic regression ...")
    p_forest = logistic_regression(x_full, coefs[...,0], coefs[...,1:])
    logging.debug("Done! Generating some plots...")
    

    return p_forest, rg_vec_subs, az_vec_subs, subs_F_r, subs_F_a

def fill_core_processing_fnf_obj(input_params_obj, pforest_equi7_fnames, mask_equi7_fnames):

    """
    Internal function called by the StackBasedProcessingFNF APP:
        
        the StackBasedProcessingFNF APP fills the structure
        to be written into the xml input file for the next APP, which is the 
        core processing for the FNF.
        The returned object "core_processing_FNF_obj" contains the paths of 
        computed data and masks in equi7 format
        
        Usage of the returned object:
            The returned object can be added to the input_params_obj and it 
            can be written to disk if needed, as (core_processing_fnf_obj is 
                                                  overwritten by this command, if already present):
                - input_params_obj.core_processing_fnf_obj = core_processing_fnf_obj
                - write_input_file(input_params_obj, input_file_xml)
    """

    core_processing_fnf_obj = core_processing_fnf(pforest_equi7_fnames, mask_equi7_fnames)

    return core_processing_fnf_obj

class CoreProcessingFNF(Task):
    """CoreProcessingFNF APP 
    
    CoreProcessingFNF APP is the second of the two sub-APPs called by 
    Forest/Non-Forest map processor.

    It performs merging ascending with descending processed stacks and
    mosaiking together differtent processed tiles of the output grid

    Attributes
    ----------
    configuration_file : str
        path of the Configuration_File.xml file

    Methods
    -------
    run( input_file_path )
        run the CoreProcessingFNF APP
    name : str
        name of the APP
        
    See Also
    --------
    biopal.dataset_query.dataset_query.dataset_query : it's the first APP to be called in the manual sequence
    StackBasedProcessingFNF : it’s the APP which prepares the input for this APP

    Examples
    -------- 
    Manual FNF chain execution
        
    >>> from biopal.dataset_query.dataset_query import dataset_query
    >>> from biopal.fnf.main_FNF import StackBasedProcessingFNF
    >>> dq_obj = dataset_query()
    >>> input_file_up1 = dq_obj.run( input_file )
    >>> sbp_obj = StackBasedProcessingFNF( config_file )
    >>> input_file_up2 = sbp_obj.run( input_file_up1 )
    >>> fnfcp_obj = CoreProcessingFNF( config_file )
    >>> fnfcp_obj.run( input_file_up2 )

    - input_file: path of the BioPAL input file
    - input_file_up1: same of input_file with also the "stack_based_processing" section   
    - input_file_up2: same of input_file_up1 with also the "core_processing_fnf" section
    - config_file: path of the BioPAL configuration file
    """

    def __init__(self, configuration_file):

        super().__init__(configuration_file)

    def _run(self, input_file):

        # FNF: Reading chains configuration files
        logging.info("FNF: Reading chains configuration files")
        check_if_path_exists(self.configuration_file, "FILE")
        conf_params_obj = parse_configuration_file(self.configuration_file)
        input_params_obj = parse_input_file(input_file)

        # managing folders
        products_folder = os.path.join(input_params_obj.output_specification.output_folder, "Products")
        temp_output_folder = os.path.join(products_folder, "temp")

        ######################## NOT STACK BASED STEPS ############################

        stacks_to_merge_dict = collect_stacks_to_be_merged(input_params_obj.stack_based_processing.stack_composition)

        try:

            logging.info("FNF: merging ascending with descending stacks....\n")

            merged_data_fnames, merging_folder = pforest_masking_and_merging(
                input_params_obj.core_processing_fnf.pforest_equi7_fnames,
                input_params_obj.core_processing_fnf.mask_equi7_fnames,
                stacks_to_merge_dict,
            )

            logging.info("...done.\n")

        except Exception as e:
            logging.error(e, exc_info=True)
            raise

        try:

            logging.info("FNF: mosaiking equi7 tiles....\n")

            output_folder = os.path.join(products_folder, "global_FNF")

            mosaiking(merging_folder, output_folder)

            logging.info("...done.\n")

        except Exception as e:
            logging.error(e, exc_info=True)
            raise

        if conf_params_obj.processing_flags.delete_temporary_files:
            shutil.rmtree(temp_output_folder)

        logging.info("FNF: Forest/NonForest map generation ended correctly.\n")

