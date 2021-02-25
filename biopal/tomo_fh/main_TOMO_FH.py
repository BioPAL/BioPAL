import os
import numpy as np
import logging
import shutil
from scipy.signal import convolve2d
from osgeo import gdal
from equi7grid.equi7grid import Equi7Grid
from equi7grid.image2equi7grid import image2equi7grid


# biomassL2 processor imports
from biopal.fh.processing_FH import heigths_masking_and_merging
from biopal.data_operations.data_operations import (
    read_and_oversample_data,
    read_and_oversample_aux_data,
    apply_dem_flattening,
    mosaiking,
)
from biopal.utility.utility_functions import (
    Task,
    choose_equi7_sampling,
    check_if_path_exists,
    check_if_geometry_auxiliaries_are_present,
    resolution_heading_correction,
    decode_unique_acquisition_id_string,
    evaluate_estimation_quality_matrix,
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
from biopal.io.data_io import tiff_formatter
from biopal.screen_calibration.screen_calibration import apply_calibration_screens
from biopal.geometry.utility_geometry import compute_and_oversample_geometry_auxiliaries


from biopal.tomo.processing_TOMO import BiomassForestHeightSKPD


class TomoForestHeight(Task):
    def __init__(
        self, configuration_file_xml, stacks_to_merge_dict, gdal_path,
    ):
        super().__init__(configuration_file_xml)
        self.stacks_to_merge_dict = stacks_to_merge_dict
        self.gdal_path = gdal_path

    def _run(self, input_file_xml):

        # Main APP #1: Stack Based Processing
        stack_based_processing_obj = StackBasedProcessingTOMOFH(self.configuration_file_xml, self.gdal_path,)

        # Run Main APP #1: Stack Based Processing
        (data_equi7_fnames, mask_equi7_fnames) = stack_based_processing_obj.run(input_file_xml)

        # Main APP #2: Core Processing
        tomo_fh_processing_obj = CoreProcessingTOMOFH(
            self.configuration_file_xml, self.stacks_to_merge_dict, data_equi7_fnames, mask_equi7_fnames,
        )

        # Run Main APP #2: AGB Core Processing
        tomo_fh_processing_obj.run(input_file_xml)


class StackBasedProcessingTOMOFH(Task):
    def __init__(
        self, configuration_file_xml, gdal_path,
    ):
        super().__init__(configuration_file_xml)
        self.gdal_path = gdal_path

    def _run(self, input_file_xml):

        ########################## INITIAL STEPS ##############################

        logging.info("TOMO FH: Reading chains configuration file")
        check_if_path_exists(self.configuration_file_xml, "FILE")
        proc_conf = parse_chains_configuration_file(self.configuration_file_xml)
        proc_inputs = parse_chains_input_file(input_file_xml)

        ### managing output folders:
        products_folder = os.path.join(proc_inputs.output_folder, "Products")
        if proc_conf.save_breakpoints:
            breakpoints_output_folder = os.path.join(products_folder, "breakpoints")
            logging.info("TOMO FH: Breakpoints will be saved into: " + breakpoints_output_folder)
            os.makedirs(breakpoints_output_folder)

        temp_output_folder = os.path.join(products_folder, "temp")
        logging.info("TOMO FH: Temporary data folder:" + temp_output_folder + "\n")
        os.makedirs(temp_output_folder)

        equi7_sampling = choose_equi7_sampling(
            proc_conf.TOMO_FH.product_resolution, proc_inputs.geographic_grid_sampling
        )
        e7g = Equi7Grid(equi7_sampling)
        logging.info("    EQUI7 Grid sampling used: {}".format(equi7_sampling))

        ########################## INITIAL STEPS END #############################

        data_equi7_fnames = {}
        mask_equi7_fnames = {}
        ########################## STACK BASED STEPS ##############################
        for unique_stack_id, acquisitions_pf_names in proc_inputs.stack_composition.items():

            # make temporary sub-directories
            temp_output_folder_gr = os.path.join(temp_output_folder, "geocoded", unique_stack_id)
            temp_output_folder_e7 = os.path.join(temp_output_folder, "equi7", unique_stack_id)
            os.makedirs(temp_output_folder_gr)
            os.makedirs(temp_output_folder_e7)

            ### load data ( and oversample if requested and if needed )
            try:
                logging.info("TOMO FH: Data loading for stack " + unique_stack_id + "; this may take a while:")

                (beta0_calibrated, master_id, raster_info, raster_info_orig,) = read_and_oversample_data(
                    proc_inputs.L1c_repository, acquisitions_pf_names, proc_conf.enable_resampling
                )
            except Exception as e:
                logging.error("TOMO FH: error during input data reading: " + str(e), exc_info=True)
                raise

            ### load or compute auxiliary data
            try:

                read_ref_h = not proc_conf.apply_calibration_screen and proc_conf.DEM_flattening
                read_cal_screens = proc_conf.apply_calibration_screen
                geometry_aux_are_present = check_if_geometry_auxiliaries_are_present(
                    proc_inputs, unique_stack_id, acquisitions_pf_names, read_ref_h=read_ref_h, read_dist=False,
                )

                if proc_conf.compute_geometry or not geometry_aux_are_present:

                    # messages for the log:
                    if proc_conf.compute_geometry:
                        logging.info("TOMO FH: calling geometry library for stack " + unique_stack_id + "\n")
                        if geometry_aux_are_present:
                            logging.warning("    geometry auxiliaries will be overwritten for stack " + unique_stack_id)
                        else:
                            logging.info("\n")
                    else:
                        logging.warning(
                            'TOMO FH: calling geometry library since AuxiliaryProductsFolder "Geometry" is empty or not complete \n'
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
                        proc_inputs.L1c_repository,
                        proc_inputs,
                        unique_stack_id,
                        acquisitions_pf_names,
                        master_id,
                        proc_conf.enable_resampling,
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

                    logging.info("TOMO FH: ...geometry auxiliaries computation done.")

                else:

                    logging.info(
                        "TOMO FH: geometry auxiliaries are provided from user, so they are now being loaded and not computed, for stack "
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
                        proc_inputs,
                        unique_stack_id,
                        acquisitions_pf_names,
                        proc_conf.enable_resampling,
                        raster_info_orig,
                        read_ref_h=read_ref_h,
                        read_dist=False,
                    )
                    logging.info("TOMO FH: ...geometry auxiliaries loading done.")

                # read the rest of auxiliaries which are notpart of the geometry library:
                if read_cal_screens:

                    logging.warning("TOMO FH: loading calibration screens \n")

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
                    logging.info("...done")

            except Exception as e:
                logging.error(
                    "TOMO FH: error during auxiliary data computation and/or loading: " + str(e), exc_info=True,
                )
                raise

            ### Screen calibration (ground steering)
            try:
                if proc_conf.apply_calibration_screen:
                    logging.info("TOMO FH: applying calibration screen...")
                    beta0_calibrated = apply_calibration_screens(
                        beta0_calibrated, raster_info, cal_screens, cal_screens_raster_info, master_id,
                    )
                    logging.info("...done.\n")

                elif proc_conf.DEM_flattening:
                    logging.info("TOM FH: DEM flattening... ")
                    beta0_calibrated = apply_dem_flattening(
                        beta0_calibrated, kz, reference_height, master_id, raster_info
                    )
                    logging.info("...done.\n")

            except Exception as e:
                logging.error(
                    "TOMO FH: error during screen calibration or DEM flattening." + str(e), exc_info=True,
                )
                raise

            ### compute mean incidence angle
            look_angle_rad = np.nanmean(off_nadir_angle_rad[master_id])  # 0.4886921905584123
            logging.info("TOMO FH: incidence angle used is {} [deg] \n".format(np.rad2deg(look_angle_rad)))

            ### mean of off nadir and slope over final resolution
            windtm_x = np.int(np.round(proc_conf.TOMO_FH.product_resolution / raster_info.pixel_spacing_az / 2) * 2 + 1)
            windtm_y = np.int(
                np.round(
                    proc_conf.TOMO_FH.product_resolution
                    / (raster_info.pixel_spacing_slant_rg / np.sin(look_angle_rad))
                    / 2
                )
                * 2
                + 1
            )

            logging.info("TOMO FH: computing mean of off nadir based on final product resolution...")

            # only mster is used in the code, do not convolve all the stacs, it is useless
            off_nadir_angle_rad[master_id] = convolve2d(
                off_nadir_angle_rad[master_id], np.ones((windtm_y, windtm_x)) / windtm_y / windtm_x, mode="same",
            )
            logging.info("...done.")

            logging.info("TOMO FH: computing mean of slope based on final product resolution...")

            slope = convolve2d(slope, np.ones((windtm_y, windtm_x)) / windtm_y / windtm_x, mode="same")
            logging.info("...done.")

            # covariance estimation window size, it may be modified by an internal flag in case of air-plane geometry
            cov_est_window_size = proc_conf.TOMO_FH.product_resolution

            if proc_conf.multilook_heading_correction:
                _, heading_deg, _, _, _, _ = decode_unique_acquisition_id_string(unique_stack_id + "_BSL_00")

                cov_est_window_size = resolution_heading_correction(cov_est_window_size, heading_deg)

            logging.info("Covariance Estimation Window Size used: {} [m]".format(cov_est_window_size))

            ### height estimation
            logging.info("TOMO FH: " + unique_stack_id + ": performing heigth estimation...")
            try:

                vertical_vector = np.arange(
                    proc_conf.vertical_range.minimum_height,
                    proc_conf.vertical_range.maximum_height + proc_conf.vertical_range.sampling,
                    proc_conf.vertical_range.sampling,
                )

                (estimated_height, power_peak, rg_vec_subs, az_vec_subs, subs_F_r, subs_F_a,) = BiomassForestHeightSKPD(
                    beta0_calibrated,
                    cov_est_window_size,
                    raster_info.pixel_spacing_slant_rg,
                    raster_info.pixel_spacing_az,
                    look_angle_rad,
                    kz,
                    vertical_vector,
                    proc_conf.TOMO_FH,
                )

                estimated_height = estimated_height * (
                    1
                    - np.tan(slope[rg_vec_subs, :][:, az_vec_subs])
                    / np.tan(off_nadir_angle_rad[master_id][rg_vec_subs, :][:, az_vec_subs])
                )

            except Exception as e:
                logging.error("FH: error during height estimation: " + str(e), exc_info=True)
                raise
            del beta0_calibrated

            ### Placemark for the quality estimation to be defined
            try:
                quality_layer_sr = evaluate_estimation_quality_matrix(estimated_height.shape)

            except Exception as e:
                logging.error("FH: error during estimation quality layer evaluation: " + str(e), exc_info=True)
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

                # geocode the FH
                logging.info(unique_stack_id + ": geocoding the estimated height...")
                data_ground = geocoding(
                    estimated_height, lon_in, lat_in, lonMeshed_out, latMeshed_out, valid_values_mask,
                )
                logging.info("...done.\n")

                # geocode the FH quality layer
                logging.info(unique_stack_id + ": geocoding the estimated height quality layer...")
                quality_layer_ground = geocoding(
                    quality_layer_sr, lon_in, lat_in, lonMeshed_out, latMeshed_out, valid_values_mask,
                )
                logging.info("...done.\n")

            except Exception as e:
                logging.error("TOMO FH: error during geocoding: " + str(e), exc_info=True)
                raise

            ### saving breakpoints
            if proc_conf.save_breakpoints:
                logging.info("TOMO FH: saving mail results (in slant range geometry) on " + breakpoints_output_folder)
                post_string = "_SR_" + unique_stack_id

                breakpoint_names = ["estimated_height" + post_string]

                save_breakpoints(breakpoints_output_folder, breakpoint_names, [estimated_height])
                logging.info("...done.\n")
            del estimated_height

            ### creating mask to exclude estimation not valid values:
            condition_curr = np.logical_and(
                data_ground > proc_conf.TOMO_FH.estimation_valid_values_limits[0],
                data_ground < proc_conf.TOMO_FH.estimation_valid_values_limits[1],
            )
            estimation_mask_ground = np.where(condition_curr, True, False)
            estimation_mask_ground[np.isnan(data_ground)] = False
            # also casted to float, for incapsulation in a geotiff
            estimation_mask_ground.astype(float)

            ### create GEOTIFF of all the layers (estimation, mask ):
            logging.info(unique_stack_id + ": formatting data to GEOTIFF...")
            try:

                data_ground_fname = os.path.join(temp_output_folder_gr, "FH" + ".tif")
                valid_values_mask_ground_fname = os.path.join(temp_output_folder_gr, "mask" + ".tif")

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

                # forest height geotiff formatting
                tiff_formatter(
                    [data_ground, quality_layer_ground],
                    data_ground_fname,
                    geotransform,
                    gdal_data_format=gdal.GDT_Float32,
                    multi_layers_tiff=True,
                )

                # forest height temporary mask geotiff formatting (this mask should be merged with the fnf mask, generated after )
                tiff_formatter(
                    estimation_mask_ground,
                    valid_values_mask_ground_fname,
                    geotransform,
                    gdal_data_format=gdal.GDT_Float32,
                )

                logging.info("...done.\n")

            except Exception as e:
                logging.error("TOMO FH: error during GEOTIFF formatting: " + str(e), exc_info=True)
                raise

            ### formatting data to EQUI7
            logging.info(unique_stack_id + ": formatting into EQUI7 grid...")
            try:

                equi7_data_outdir = os.path.join(temp_output_folder_e7, "data")
                equi7_mask_outdir = os.path.join(temp_output_folder_e7, "mask")

                # in general from here the Equi7 can output multiple tiles, which file names are stored in the output list ( wrapped here in a dict for the stack )
                data_equi7_fnames[unique_stack_id] = image2equi7grid(
                    e7g,
                    data_ground_fname,
                    equi7_data_outdir,
                    gdal_path=self.gdal_path,
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
                    gdal_path=self.gdal_path,
                    inband=None,
                    subgrid_ids=None,
                    accurate_boundary=False,
                    withtilenamesuffix=False,
                    resampling_type="bilinear",
                    tile_nodata=np.float(0),
                )

            except Exception as e:
                logging.error("TOMO FH: error during EQUI7 formatting: " + str(e), exc_info=True)
                raise

            if not os.path.exists(data_equi7_fnames[unique_stack_id][0]):
                error_message = "EQUI7 grid has not been generated, output is absent "
                logging.error(error_message)
                raise

            logging.info("...done.\n")
        ######################## STACK BASED STEPS END. ###########################
        return (
            data_equi7_fnames,
            mask_equi7_fnames,
        )


class CoreProcessingTOMOFH(Task):
    def __init__(
        self, configuration_file_xml, stacks_to_merge_dict, data_equi7_fnames, mask_equi7_fnames,
    ):
        super().__init__(configuration_file_xml)
        self.stacks_to_merge_dict = stacks_to_merge_dict
        self.data_equi7_fnames = data_equi7_fnames
        self.mask_equi7_fnames = mask_equi7_fnames

    def _run(self, input_file_xml):
        ######################## NOT STACK BASED STEPS ############################

        ### managing output folders:
        check_if_path_exists(self.configuration_file_xml, "FILE")
        proc_conf = parse_chains_configuration_file(self.configuration_file_xml)
        proc_inputs = parse_chains_input_file(input_file_xml)
        products_folder = os.path.join(proc_inputs.output_folder, "Products")
        temp_output_folder = os.path.join(products_folder, "temp")

        try:

            logging.info("TOMO FH: merging ascending with descending stacks....\n")

            merged_data_fnames, merging_folder = heigths_masking_and_merging(
                self.data_equi7_fnames, self.mask_equi7_fnames, self.stacks_to_merge_dict
            )

            logging.info("...done.\n")

        except Exception as e:
            logging.error(e, exc_info=True)
            raise

        try:

            logging.info("TOMO FH: mosaiking equi7 tiles....\n")

            output_folder = os.path.join(products_folder, "global_FH")

            mosaiking(merging_folder, output_folder)

            logging.info("...done.\n")

        except Exception as e:
            logging.error(e, exc_info=True)
            raise

        if proc_conf.delete_temporary_files:
            shutil.rmtree(temp_output_folder)

        logging.info("TOMO FH: Forest Height estimation ended correctly.\n")
