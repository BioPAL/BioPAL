import os
import logging
import pyproj
import shutil
import numpy as np
from osgeo import gdal
from gdalconst import GA_ReadOnly
from scipy.interpolate import interp1d
from skimage.filters.rank import majority as majority_filter
from equi7grid.equi7grid import Equi7Grid
from equi7grid.image2equi7grid import image2equi7grid
from biopal.io.data_io import (
    read_data,
    read_ecef_grid,
    read_auxiliary_multi_channels,
    read_auxiliary_single_channel,
    tandemx_fnf_read,
    tiff_formatter,
    readBiomassHeader,
)
from biopal.io.xml_io import raster_info
from biopal.utility.constants import OVERSAMPLING_FACTOR
from biopal.utility.utility_functions import get_equi7_fnf_tiff_names
from arepytools.io.productfolder import ProductFolder


def data_oversample(data, oversampling_factor, raster_info_obj):

    rg_ratio = np.floor(raster_info_obj.resolution_m_slant_rg / raster_info_obj.pixel_spacing_slant_rg)
    az_ratio = np.floor(raster_info_obj.resolution_m_az / raster_info_obj.pixel_spacing_az)

    rg_oversampling_flag = False
    az_oversampling_flag = False

    if rg_ratio < 2:
        rg_oversampling_flag = True
        num_samples_out = raster_info_obj.num_samples * oversampling_factor
    else:
        num_samples_out = raster_info_obj.num_samples
        pixel_spacing_slant_rg_out = raster_info_obj.pixel_spacing_slant_rg

    if az_ratio < 2:
        az_oversampling_flag = True
        num_lines_out = raster_info_obj.num_lines * oversampling_factor
    else:
        num_lines_out = raster_info_obj.num_lines
        pixel_spacing_az_out = raster_info_obj.pixel_spacing_az

    if rg_oversampling_flag or az_oversampling_flag:
        logging.info("    Oversampling needed:")

        if not type(data) is dict:
            # slope and reference_heights cases  (input data is value matrix)

            if rg_oversampling_flag:

                logging.info("        range   oversampling of auxiliary data...")

                data, pixel_spacing_slant_rg_out = data_oversample_core(
                    data, 0, raster_info_obj.pixel_spacing_slant_rg, oversampling_factor
                )

            if az_oversampling_flag:

                logging.info("        azimuth   oversampling of auxiliary data...")

                data, pixel_spacing_az_out = data_oversample_core(
                    data, 1, raster_info_obj.pixel_spacing_az, oversampling_factor
                )

        else:

            for data_key, data_extracted in data.items():

                if not type(data_extracted) is dict:

                    # KZ, off_nadir and ECEFGRID case (input data is a dict of values)
                    if rg_oversampling_flag:

                        logging.info("        range   oversampling of " + data_key + "...")

                        data[data_key], pixel_spacing_slant_rg_out = data_oversample_core(
                            data_extracted, 0, raster_info_obj.pixel_spacing_slant_rg, oversampling_factor,
                        )

                    if az_oversampling_flag:
                        logging.info("        azimuth oversampling of " + data_key + "...")

                        data_extracted = data[data_key]
                        data[data_key], pixel_spacing_az_out = data_oversample_core(
                            data_extracted, 1, raster_info_obj.pixel_spacing_az, oversampling_factor
                        )

                else:

                    # Beta0 case (input data is a dict of dict with values)
                    for pol_key, data_pol in data_extracted.items():

                        if rg_oversampling_flag:

                            logging.info(
                                "        range   oversampling of " + data_key + " , polarization " + pol_key + "..."
                            )

                            (data[data_key][pol_key], pixel_spacing_slant_rg_out,) = data_oversample_core(
                                data_pol, 0, raster_info_obj.pixel_spacing_slant_rg, oversampling_factor,
                            )

                        if az_oversampling_flag:
                            logging.info(
                                "        azimuth oversampling of " + data_key + " , polarization " + pol_key + "..."
                            )

                            data_pol = data[data_key][pol_key]
                            data[data_key][pol_key], pixel_spacing_az_out = data_oversample_core(
                                data_pol, 1, raster_info_obj.pixel_spacing_az, oversampling_factor
                            )

        logging.info("    oversampling done.\n")

    return data, num_samples_out, pixel_spacing_slant_rg_out, num_lines_out, pixel_spacing_az_out


def data_oversample_core(data, axis_index, pixel_spacing_in, oversampling_factor):

    # original size of the axis to be interpolated
    axis_len_in = data.shape[axis_index]

    # pixel spacing after the interpolation
    pixel_spacing_out = pixel_spacing_in / oversampling_factor

    # input original axis
    max_val_in = axis_len_in * pixel_spacing_in
    ax_in = np.arange(0, max_val_in, pixel_spacing_in)

    # output interpolated axis
    max_val_out = axis_len_in * oversampling_factor * pixel_spacing_out
    ax_out = np.arange(0, max_val_out, pixel_spacing_out)

    # interpolation of data along axis_index
    interp_fun = interp1d(ax_in, data, axis=axis_index, bounds_error=0)
    data = interp_fun(ax_out)

    return data, pixel_spacing_out


def read_and_oversample_data(L1c_repository, acquisitions_pf_names, enable_resampling):
    # this function calls the read_data followed by the data_oversample
    # (which oversamples only when needed and if enabled )

    beta0_calibrated = {}
    for pf_name in acquisitions_pf_names:
        logging.info("    loading " + pf_name + "...")
        # every data in the current stack has same spacing, same resolution and same number of samples and lines
        (
            beta0_calibrated[pf_name],
            num_samples,
            num_lines,
            pixel_spacing_slant_rg,
            pixel_spacing_az,
            carrier_frequency_hz,
            range_bandwidth_hz,
            master_id,
            lines_start_utc,
        ) = read_data(L1c_repository, pf_name)

    ### ALL chains: oversampling data:
    # needed whenever computing covariance or detected data (as notch for AGB)
    (_, _, _, _, _, _, resolution_m_slant_rg, resolution_m_az, sensor_velocity,) = readBiomassHeader(
        ProductFolder(os.path.join(L1c_repository, acquisitions_pf_names[0]), "r"), 0
    )

    # backup needed for the resampling in auxiliary data
    raster_info_orig = raster_info(
        num_samples,
        num_lines,
        pixel_spacing_slant_rg,
        pixel_spacing_az,
        resolution_m_slant_rg,
        resolution_m_az,
        carrier_frequency_hz,
        range_bandwidth_hz,
        lines_start_utc,
    )

    if enable_resampling:

        (beta0_calibrated, num_samples, pixel_spacing_slant_rg, num_lines, pixel_spacing_az,) = data_oversample(
            beta0_calibrated, OVERSAMPLING_FACTOR, raster_info_orig,
        )

        logging.info("all data loaded.\n")

    # filling output structure:
    raster_info_os = raster_info(
        num_samples,
        num_lines,
        pixel_spacing_slant_rg,
        pixel_spacing_az,
        resolution_m_slant_rg,
        resolution_m_az,
        carrier_frequency_hz,
        range_bandwidth_hz,
        lines_start_utc,
    )

    return beta0_calibrated, master_id, raster_info_os, raster_info_orig


def read_and_oversample_aux_data(
    file_names,
    stack_id,
    acquisitions_pf_names,
    enable_resampling,
    raster_info,
    read_ecef=True,
    read_off_nadir=True,
    read_slope=True,
    read_kz=True,
    read_ref_h=True,
    read_dist=True,
    read_average_cov=False,
    read_cal_screens=False,
    read_FH=False,
    read_reference_agb=False,
    read_sys_dec_fun=False,
):

    # this function calls the auxilary data readers followed by the data_oversample
    # (which oversamples only when needed and if enabled )
    # it reads only the auxiliary in the current stack_id

    # for TOMO we need another auxiliary: phases to be extreacted

    ecef_grid = None
    kz = None
    off_nadir_angle_rad = None
    reference_height = None
    R = None
    slope = None
    average_covariance = None
    calibration_screens = None
    cal_screens_raster_info = None
    forest_height = None
    reference_agb = None
    system_decorr_fun = None

    if read_ecef:
        logging.info("Loading auxiliary data: ECEFGRID...")
        pf_name = os.path.basename(file_names.ECEF_grid_file_names[stack_id])
        folder_name = os.path.dirname(file_names.ECEF_grid_file_names[stack_id])
        ecef_grid = read_ecef_grid(folder_name, pf_name)

        if not ecef_grid is None and enable_resampling:
            ecef_grid = data_oversample(ecef_grid, OVERSAMPLING_FACTOR, raster_info)[0]
        if not ecef_grid is None:
            logging.info("...ECEFGRID read.\n")
        else:
            logging.warning(
                'Since "ECEFGRID/'
                + stack_id
                + '" folder is missing into AuxiliaryProductsFolder, geometric library needs to be called.\n'
            )

    if read_kz:
        logging.info("Loading auxiliary data: KZ...")
        pf_name = os.path.basename(file_names.kz_file_names[stack_id])
        folder_name = os.path.dirname(file_names.kz_file_names[stack_id])
        kz = read_auxiliary_multi_channels(folder_name, pf_name, acquisitions_pf_names)

        if not kz is None and enable_resampling:
            kz = data_oversample(kz, OVERSAMPLING_FACTOR, raster_info)[0]
        if not kz is None:
            logging.info("...KZ read.\n")
        else:
            logging.warning(
                'Since "KZ/'
                + stack_id
                + '" folder is missing into AuxiliaryProductsFolder, geometric library needs to be called.\n'
            )

    if read_off_nadir:
        logging.info("Loading auxiliary data: off nadir angles...")
        pf_name = os.path.basename(file_names.off_nadir_angle_file_names[stack_id])
        folder_name = os.path.dirname(file_names.off_nadir_angle_file_names[stack_id])
        off_nadir_angle_rad = read_auxiliary_multi_channels(folder_name, pf_name)

        if not off_nadir_angle_rad is None and enable_resampling:
            off_nadir_angle_rad = data_oversample(off_nadir_angle_rad, OVERSAMPLING_FACTOR, raster_info)[0]
        if not off_nadir_angle_rad is None:
            logging.info("...off nadir angles read.\n")
        else:
            logging.warning(
                'Since "OffNadirAngles/'
                + stack_id
                + '" folder is missing into AuxiliaryProductsFolder, geometric library needs to be called.\n'
            )

    if read_ref_h:
        logging.info(
            "Loading auxiliary data: reference height..."
        )  # not needed, use just for comparison with the estimation
        pf_name = os.path.basename(file_names.reference_height_file_names[stack_id])
        folder_name = os.path.dirname(file_names.reference_height_file_names[stack_id])
        reference_height = read_auxiliary_single_channel(
            folder_name, pf_name
        )  # it is the dtm in slant_range-azimuth reference

        if not reference_height is None and enable_resampling:
            reference_height = data_oversample(reference_height, OVERSAMPLING_FACTOR, raster_info)[0]
        if not reference_height is None:
            logging.info("...reference height read.\n")
        else:
            logging.warning(
                'Since "ReferenceHeight/'
                + stack_id
                + '" folder is missing into AuxiliaryProductsFolder, geometric library needs to be called.\n'
            )

    if read_dist:
        logging.info("Loading auxiliary data: slant range distances...")
        pf_name = os.path.basename(file_names.slant_range_distances_file_names[stack_id])
        folder_name = os.path.dirname(file_names.slant_range_distances_file_names[stack_id])
        R = read_auxiliary_multi_channels(folder_name, pf_name)

        if not R is None and enable_resampling:
            R = data_oversample(R, OVERSAMPLING_FACTOR, raster_info)[0]
        if not R is None:
            logging.info("...slant range distances read.\n")
        else:
            logging.warning(
                'Since "SlantRangeDistances/'
                + stack_id
                + '" folder is missing into AuxiliaryProductsFolder, geometric library needs to be called.\n'
            )

    if read_slope:
        logging.info("Loading auxiliary data: slope...")
        pf_name = os.path.basename(file_names.slope_file_names[stack_id])
        folder_name = os.path.dirname(file_names.slope_file_names[stack_id])
        slope = read_auxiliary_single_channel(folder_name, pf_name)

        if not slope.any() is None and enable_resampling:
            slope = data_oversample(slope, OVERSAMPLING_FACTOR, raster_info)[0]
        if not slope is None:
            logging.info("...slope read.\n")
        else:
            logging.warning(
                'Since "Slopes/'
                + stack_id
                + '" folder is missing into AuxiliaryProductsFolder, geometric library needs to be called.\n'
            )

    if read_average_cov:
        logging.info("Loading auxiliary data: Average covariance matrix...")
        pf_name = os.path.basename(file_names.average_covariance_folder)
        folder = os.path.dirname(file_names.average_covariance_folder)
        average_covariance = read_auxiliary_multi_channels(folder, pf_name, acquisitions_pf_names)

        if not average_covariance is None and enable_resampling:
            average_covariance = data_oversample(average_covariance, OVERSAMPLING_FACTOR, raster_info)[0]
        if not average_covariance is None:
            logging.info("... Average covariance matrix read.\n")

    if read_cal_screens:
        logging.info("Loading auxiliary data: Calibration Screens...")
        pf_name = os.path.basename(file_names.calibration_screens_file_names[stack_id])
        folder = os.path.dirname(file_names.calibration_screens_file_names[stack_id])
        calibration_screens, cal_screens_raster_info = read_auxiliary_multi_channels(
            folder, pf_name, acquisitions_pf_names, read_raster_info=True
        )

        if not calibration_screens is None and enable_resampling:
            calibration_screens = data_oversample(calibration_screens, OVERSAMPLING_FACTOR, raster_info)[0]
        if not calibration_screens is None:
            logging.info("... Calibration Screens read.\n")

    if read_FH:
        logging.info("Loading auxiliary data: forest estimated height..")
        pf_name = os.path.basename(file_names.forest_height_folder)
        folder = os.path.dirname(file_names.forest_height_folder)
        forest_height = read_auxiliary_single_channel(folder, pf_name, acquisitions_pf_names)

        if not forest_height is None and enable_resampling:
            forest_height = data_oversample(forest_height, OVERSAMPLING_FACTOR, raster_info)[0]
        logging.info("...done.\n")
        if not forest_height is None:
            logging.info("...forest estimated height read.\n")

    if read_reference_agb:
        logging.info("Loading auxiliary data: reference agb..")
        pf_name = os.path.basename(file_names.reference_agb_folder)
        folder = os.path.dirname(file_names.reference_agb_folder)
        reference_agb = read_auxiliary_single_channel(folder, pf_name, acquisitions_pf_names)

        if not reference_agb is None and enable_resampling:
            reference_agb = data_oversample(reference_agb, OVERSAMPLING_FACTOR, raster_info)[0]
        if not reference_agb is None:
            logging.info("...reference agb read.\n")

    if read_sys_dec_fun:
        logging.info("Loading auxiliary data: system decorrelation function...")
        pf_name = os.path.basename(file_names.system_decorrelation_fun_folder)
        folder = os.path.dirname(file_names.system_decorrelation_fun_folder)
        system_decorr_fun = read_auxiliary_single_channel(folder, pf_name, acquisitions_pf_names)

        if not system_decorr_fun is None and enable_resampling:
            system_decorr_fun = data_oversample(system_decorr_fun, OVERSAMPLING_FACTOR, raster_info)[0]
        if not system_decorr_fun is None:
            logging.info("...system decorrelation function read.\n")

    return (
        ecef_grid,
        off_nadir_angle_rad,
        slope,
        kz,
        reference_height,
        R,
        average_covariance,
        calibration_screens,
        cal_screens_raster_info,
        forest_height,
        reference_agb,
        system_decorr_fun,
    )


def fnf_tandemx_load_filter_equi7format(
    forest_mask_catalogue_folder,
    e7g,
    product_resolution,
    output_folder,
    gdal_path,
    geographic_boundaries,
    time_tag_mjd_initial,
    dummy_mask_flag=False,
):
    # this function loads from the fnf catalog folder, the times which matches the boundaries
    # the output is also casted to float to be ready for incapsulation within a geotiff

    fnf_tile_loaded_list, time_tag_mjd_tandemx, geotransform_list = tandemx_fnf_read(
        forest_mask_catalogue_folder, geographic_boundaries
    )

    # majority filter
    fnf_mask_list = []
    for idx_fnf, fnf in enumerate(fnf_tile_loaded_list):

        len_x, len_y = fnf.shape
        # two longitide and latitude points in middle of the fnf mask:

        # geotransform = [ lon_start, lon_step, 0, lat_start, 0, lat_step]
        lon0 = geotransform_list[idx_fnf][0] + geotransform_list[idx_fnf][1] * len_x / 2
        lon1 = lon0 + geotransform_list[idx_fnf][1]
        lat0 = geotransform_list[idx_fnf][3] + geotransform_list[idx_fnf][5] * len_x / 2
        lat1 = lat0 + geotransform_list[idx_fnf][5]

        geod = pyproj.Geod(ellps="WGS84")

        _, _, fnf_spacing_lon = geod.inv(lon0, lat0, lon1, lat0)
        _, _, fnf_spacing_lat = geod.inv(lon0, lat0, lon0, lat1)

        windtm_x = np.int(np.round(product_resolution / fnf_spacing_lon / 2) * 2 + 1)
        windtm_y = np.int(np.round(product_resolution / fnf_spacing_lat / 2) * 2 + 1)

        selem = np.ones((windtm_x, windtm_y), dtype=int)
        if not dummy_mask_flag:
            fnf_filtered = majority_filter(fnf.astype(int), selem)

            fnf_filtered = fnf_filtered.astype(fnf.dtype)
        else:
            fnf_filtered = np.ones(fnf.shape, dtype=fnf.dtype)

        # round the mask to quantize it because it is boolean: NOTE, do not cast to bool  otherwise the NaN will become True!
        fnf_filtered = np.round(fnf_filtered)

        fnf_mask_list.append(fnf_filtered)

    if time_tag_mjd_tandemx > time_tag_mjd_initial:
        error_msg = 'Input FNF Mask (TANDEM-X) date of "{}" is bigger than input stack data minimum date of "{}"'.format(
            str(time_tag_mjd_tandemx), str(time_tag_mjd_initial)
        )
        logging.error(error_msg)
        raise ValueError(error_msg)

    # conversion step 2: fnf mask geotiff formatting: one single layer geotiff for each fnf "tamdem-x tile" (they will be merged togheter after equi7 conversion )
    fnf_mask_ground_dir_name = os.path.join(output_folder, "fnf_ground")
    fnf_mask_ground_dir_names = tiff_formatter(
        fnf_mask_list,
        fnf_mask_ground_dir_name,
        geotransform_list,
        gdal_data_format=gdal.GDT_Float32,
        multi_layers_tiff=False,
        time_tag=str(time_tag_mjd_tandemx),
    )

    # conversion step 3: special case for the equi7 convertion of all the fnf masks found
    equi7_fnf_mask_tempdir_tiles_not_merged = (
        []
    )  # thos will contain all the fnf names in equi7: all fnf tiles, all equi7 tiles
    for idx, fnf_name_curr in enumerate(fnf_mask_ground_dir_names):
        equi7_fnf_mask_parent_tempdir = os.path.join(output_folder, "input_fnf_not_merged_T{}".format(idx))
        equi7_out_name_curr = image2equi7grid(
            e7g,
            fnf_name_curr,
            equi7_fnf_mask_parent_tempdir,
            gdal_path=gdal_path,
            inband=None,
            subgrid_ids=None,
            accurate_boundary=False,
            withtilenamesuffix=False,
            resampling_type="bilinear",
            tile_nodata=np.float(0),
        )
        equi7_fnf_mask_tempdir_tiles_not_merged.extend(equi7_out_name_curr)

    # conversion step 4: merging fnf tiles
    equi7_fnf_mask_parent_tempdir = os.path.join(output_folder, "input_fnf")
    # managing final mask creation:
    equi7_fnf_mask_fnames = fnf_equi7_masks_merging(
        equi7_fnf_mask_tempdir_tiles_not_merged, equi7_fnf_mask_parent_tempdir
    )

    logging.info("...fnf mask conversion into equi7 done.\n")

    return equi7_fnf_mask_fnames


def fnf_equi7_load_filter_equi7format(
    forest_mask_catalogue_folder, e7g_curr_chain, product_resolution, output_folder, gdal_path
):

    # get the names of input equi7 fnf masks
    equi7_fnf_mask_tiff_names_list, equi7_subgrid_name = get_equi7_fnf_tiff_names(forest_mask_catalogue_folder)

    subgrid_code = equi7_subgrid_name[6:8]

    if not equi7_fnf_mask_tiff_names_list:
        error_str = "Cannot find any Equi7 forest mask in input"
        logging.error(error_str)
        raise ValueError(error_str)

    # temporary folder for the input mask conversion
    equi7_fnf_mask_parent_tempdir = os.path.join(output_folder, "input_fnf")

    # output folder of the mask that will be used in the current chainflow:
    equi7_fnf_mask_parent_dir = os.path.join(equi7_fnf_mask_parent_tempdir, "equi7")

    if not os.path.exists(equi7_fnf_mask_parent_dir):
        os.makedirs(equi7_fnf_mask_parent_dir)
    if not os.path.exists(equi7_fnf_mask_parent_tempdir):
        os.makedirs(equi7_fnf_mask_parent_tempdir)

    equi7_fnf_mask_tiff_names_list_out = []

    # for each input Equi7 tile, load, filter and save with current equi7 grid sampling
    for idx, equi7_fnf_mask_tiff_name in enumerate(equi7_fnf_mask_tiff_names_list):

        # 1) filtering
        fnf_driver = gdal.Open(equi7_fnf_mask_tiff_name, GA_ReadOnly)
        fnf_data = fnf_driver.ReadAsArray()
        geotransform = fnf_driver.GetGeoTransform()
        projection_equi7 = fnf_driver.GetProjection()
        len_x, len_y = fnf_data.shape
        east_axis = geotransform[0] + geotransform[1] * np.arange(fnf_driver.RasterXSize)
        north_axis = geotransform[3] + geotransform[5] * np.arange(fnf_driver.RasterYSize)
        fnf_driver = None

        windtm_x = np.int(np.round(product_resolution / geotransform[1] / 2) * 2 + 1)
        windtm_y = np.int(np.round(product_resolution / abs(geotransform[5]) / 2) * 2 + 1)

        selem = np.ones((windtm_x, windtm_y), dtype=int)
        fnf_filtered = majority_filter(fnf_data.astype(int), selem)
        len_x, len_y = fnf_filtered.shape
        fnf_filtered = fnf_filtered.astype(fnf_data.dtype)
        # round the mask to quantize it because it is boolean: NOTE, do not cast to bool  otherwise the NaN will become True!
        fnf_filtered = np.round(fnf_filtered)

        # 2) intermediate saving to tiff
        intermediate_filtered_tiff_name = os.path.join(equi7_fnf_mask_parent_tempdir, "fnf_T{}.tif".format(idx))

        intermediate_filtered_tiff_name = tiff_formatter(
            fnf_filtered,
            intermediate_filtered_tiff_name,
            geotransform,
            gdal_data_format=gdal.GDT_Float32,
            projection=projection_equi7,
        )

        # create the equi7grid object with the input mask sampling
        equi7_sampling_mask_in = int(geotransform[1])
        try:
            e7g_in = Equi7Grid(equi7_sampling_mask_in)
        except Exception as e:
            logging.error(
                "Input FNF Equi7 mask geotransform is not in the correct format :{} ".format(subgrid_code) + str(e),
                exc_info=True,
            )
            raise

        # [(left, lower), (right, upper)]
        try:
            lon_min, lat_min = getattr(e7g_in, subgrid_code).xy2lonlat(min(east_axis), min(north_axis))
            lon_max, lat_max = getattr(e7g_in, subgrid_code).xy2lonlat(max(east_axis), max(north_axis))
        except Exception as e:
            logging.error(
                'Cannot recognize input FNF Equi7 mask "{}" sub-grid folder name :'.format(subgrid_code) + str(e),
                exc_info=True,
            )
            raise

        bbox = [(lon_min, lat_min), (lon_max, lat_max)]
        ftiles = e7g_curr_chain.search_tiles_in_roi(bbox=bbox)

        out_equi7_folder_done = image2equi7grid(
            e7g_curr_chain,
            intermediate_filtered_tiff_name,
            equi7_fnf_mask_parent_dir,
            gdal_path=gdal_path,
            ftiles=ftiles,
            accurate_boundary=False,
            withtilenamesuffix=False,
            tile_nodata=np.nan,
        )

        for idx, equi7_fnf_mask_tiff_name in enumerate(out_equi7_folder_done):

            fnf_driver = gdal.Open(equi7_fnf_mask_tiff_name, GA_ReadOnly)
            fnf_data = fnf_driver.ReadAsArray()
            fnf_driver = None

            if np.sum(np.isnan(fnf_data)) == fnf_data.shape[0] * fnf_data.shape[1]:
                shutil.rmtree(os.path.dirname(equi7_fnf_mask_tiff_name))
            else:
                equi7_fnf_mask_tiff_names_list_out.append(equi7_fnf_mask_tiff_name)

    if not equi7_fnf_mask_tiff_names_list_out:
        error_str = "Cannot find a valid input FNF Equi7 mask"
        logging.error(error_str)
        raise ImportError(error_str)

    return equi7_fnf_mask_tiff_names_list_out


def fnf_equi7_masks_merging(fnf_mask_equi7_fnames_whole, temp_output_folder):
    # there are many fnf mask tiles
    # for each fnf mask tile there are many Equi7 tiles
    # equi7 tiles may be grouped in different sub grids

    # retrive the name of tge "EQUI7 sub grid (there will be just one)
    equi7_subgrid_name = os.path.basename(os.path.dirname(os.path.dirname(fnf_mask_equi7_fnames_whole[0])))
    # also prepare the empty dictionary for equi7 tile names
    equi7_tile_names = []

    # retrive the names of all the used "EQUI7 tiles" and store in the equi7_tile_names: take the names from the generated EQUI7 folders for the fnf_mask_equi7
    for (
        fnf_tile_name
    ) in (
        fnf_mask_equi7_fnames_whole
    ):  # fnf_mask_equi7_fnames is a list of lists: fnf_mask_equi7_fnames[ list of fnf tiles ][ list of equi7 tiles ]
        curr_equi7_tile_name = os.path.basename(os.path.dirname(fnf_tile_name))
        if not curr_equi7_tile_name in equi7_tile_names:
            equi7_tile_names.append(curr_equi7_tile_name)

    merged_fnf_equi7_names = []
    # we are in the current EQUI7 sub-grid, cycle over all the
    for equi7_tile_name in equi7_tile_names:
        # search in all the fnf folders for same names
        out_fnf_equi7_name = os.path.join(
            temp_output_folder, "equi7_fnf_tiles_merged", equi7_subgrid_name, equi7_tile_name, "fnf.tif",
        )

        # search all the equi7 tile name in each fnf tile and sum togheter
        for fnf_tiff_full_path_curr in fnf_mask_equi7_fnames_whole:
            if equi7_tile_name in fnf_tiff_full_path_curr:

                input_image_driver = gdal.Open(fnf_tiff_full_path_curr, GA_ReadOnly)
                loaded_equi7_tile = input_image_driver.ReadAsArray()

                geotransform = input_image_driver.GetGeoTransform()
                projection = input_image_driver.GetProjection()

                # replace:
                loaded_equi7_tile = np.round(loaded_equi7_tile)
                loaded_equi7_tile[loaded_equi7_tile == 2] = 0
                loaded_equi7_tile[loaded_equi7_tile == 3] = 0

                # mean:
                if not "data_sum" in locals():
                    data_sum = loaded_equi7_tile
                else:
                    where_loaded_only_is_not_zero = np.logical_and(data_sum == 0, loaded_equi7_tile != 0)
                    where_both_are_not_zero = np.logical_and(data_sum != 0, loaded_equi7_tile != 0)
                    data_sum[where_loaded_only_is_not_zero] = loaded_equi7_tile[where_loaded_only_is_not_zero]
                    data_sum[where_both_are_not_zero] = (
                        loaded_equi7_tile[where_both_are_not_zero] + data_sum[where_both_are_not_zero]
                    ) / 2

        # round:
        data_sum.astype(float)

        # write the tiff as an equi7
        Nx, Ny = data_sum.shape

        driver = gdal.GetDriverByName("GTiff")
        os.makedirs(os.path.dirname(out_fnf_equi7_name))
        outdata = driver.Create(out_fnf_equi7_name, Ny, Nx, 1, gdal.GDT_Float32)
        outdata.SetGeoTransform(geotransform)
        outdata.SetProjection(projection)
        outdata.GetRasterBand(1).WriteArray(data_sum)
        outdata.FlushCache()  # saves to disk
        outdata = None

        del data_sum

        merged_fnf_equi7_names.append(out_fnf_equi7_name)

    return merged_fnf_equi7_names


def fnf_and_validity_masks_merging(merged_fnf_equi7_names, validity_mask_equi7_fnames, temp_output_folder, stack_id):

    final_mask_equi7_names = []

    # retrive the name of tge "EQUI7 sub grid (there will be just one)
    equi7_subgrid_name = os.path.basename(os.path.dirname(os.path.dirname(merged_fnf_equi7_names[0])))
    # also prepare the empty dictionary for equi7 tile names
    equi7_tile_names = []

    # retrive the names of all the used "EQUI7 tiles" and store in the equi7_tile_names: take the names from the generated EQUI7 folders for the fnf_mask_equi7
    for (
        fnf_equi7_name
    ) in (
        merged_fnf_equi7_names
    ):  # fnf_mask_equi7_fnames is a list of lists: fnf_mask_equi7_fnames[ list of fnf tiles ][ list of equi7 tiles ]
        curr_equi7_tile_name = os.path.basename(os.path.dirname(fnf_equi7_name))
        if not curr_equi7_tile_name in equi7_tile_names:
            equi7_tile_names.append(curr_equi7_tile_name)

    # search all the equi7 tile name in each fnf tile and sum togheter
    for equi7_tile_name in equi7_tile_names:
        fnf_mask_name = [name for name in merged_fnf_equi7_names if equi7_tile_name in name]
        validity_mask_name = [name for name in validity_mask_equi7_fnames if equi7_tile_name in name]
        # both should exist:
        if fnf_mask_name and validity_mask_name:

            fnf_mask_name = fnf_mask_name[0]
            validity_mask_name = validity_mask_name[0]
            out_final_equi7_name = os.path.join(
                temp_output_folder, "equi7_global_fnf", stack_id, equi7_subgrid_name, equi7_tile_name, "fnf.tif",
            )

            input_image_driver_fnf = gdal.Open(fnf_mask_name, GA_ReadOnly)
            loaded_equi7_fnf = input_image_driver_fnf.ReadAsArray()

            input_image_driver_val = gdal.Open(validity_mask_name, GA_ReadOnly)
            loaded_equi7_val = input_image_driver_val.ReadAsArray()

            final_mask = np.logical_and(loaded_equi7_fnf, loaded_equi7_val)

            final_mask.astype(int)

            # write the tiff as an equi7
            Nx, Ny = final_mask.shape

            geotransform = input_image_driver_fnf.GetGeoTransform()
            projection = input_image_driver_fnf.GetProjection()

            driver = gdal.GetDriverByName("GTiff")
            os.makedirs(os.path.dirname(out_final_equi7_name))
            outdata = driver.Create(out_final_equi7_name, Ny, Nx, 1, gdal.GDT_Float32)
            outdata.SetGeoTransform(geotransform)
            outdata.SetProjection(projection)
            outdata.GetRasterBand(1).WriteArray(final_mask)
            outdata.FlushCache()  ##saves to disk!!
            outdata = None

            final_mask_equi7_names.append(out_final_equi7_name)

    return final_mask_equi7_names


def mosaiking(equi7_main_full_folder, mosaiking_out_folder):
    # equi7_main_full_folder is the folder containing the temporary equi7 generated by the processor:
    # the final mosaiking will prodice the output folder, containing only Equi7 tiles, with no mention to stack IDs
    # this function uses listdir to retrive files and folder names from equi7_main_full_folder

    # procedure is the following:
    # first iteration is to retrive all the EQUI7 TILE names (as "E045N050T1" )
    # second iterates foe each of the above EQUI7 TILE name and does:
    #    inner cycle over all the stack folders to get the names of the tiff which matches the EQUI7 tile name
    #

    # each input in data_equi7_fnames is a tiff with two layers
    # first later is data
    data_layer_index = 1
    # second layer is quality
    quality_layer_index = 2

    # get the names of the stack or merged-stacks folder generated by the processor
    equi7_stack_folder_names = os.listdir(equi7_main_full_folder)

    # retrive the name of tge "EQUI7 sub grid (there will be just one)
    equi7_subgrid_folder_name = os.listdir(os.path.join(equi7_main_full_folder, equi7_stack_folder_names[0]))[0]

    ### first iteration: retrive all the EQUI7 TILE names (as "E045N050T1" )
    equi7_tile_names = []
    # retrive the names of all the used "EQUI7 tiles" and store in the equi7_tile_names
    for equi7_stack_folder_name in equi7_stack_folder_names:

        equi7_subgrid_full_folder = os.path.join(
            equi7_main_full_folder, equi7_stack_folder_name, equi7_subgrid_folder_name
        )

        equi7_tile_folder_names = os.listdir(equi7_subgrid_full_folder)
        for equi7_tile_name in equi7_tile_folder_names:
            if not equi7_tile_name in equi7_tile_names:
                equi7_tile_names.append(equi7_tile_name)

    ### second iteration
    out_mosaiked_equi7_tiff_names = []

    # cycle to write a modaiked output for each Equi7 tile name:
    for equi7_tile_name in equi7_tile_names:

        # inner cycle to get all the tiff names of the same equi7 ile:
        equi7_tile_tiff_full_names = []
        for equi7_stack_folder_name in equi7_stack_folder_names:

            equi7_subgrid_folder = os.path.join(
                equi7_main_full_folder, equi7_stack_folder_name, equi7_subgrid_folder_name
            )
            equi7_tile_folder_names = os.listdir(equi7_subgrid_folder)

            for equi7_tile_name_curr in equi7_tile_folder_names:

                if equi7_tile_name_curr == equi7_tile_name:
                    equi7_tile_tiff_name = os.listdir(os.path.join(equi7_subgrid_folder, equi7_tile_name))[0]
                    equi7_tile_tiff_full_names.append(
                        os.path.join(equi7_subgrid_folder, equi7_tile_name, equi7_tile_tiff_name)
                    )

        # now cycle over all the tiff data corresponding to same equi7 tile, amd performthe mosaicking (mean them together)
        out_tile_equi7_full_folder = os.path.join(mosaiking_out_folder, equi7_subgrid_folder_name, equi7_tile_name)
        for equi7_tiff_name in equi7_tile_tiff_full_names:

            # read current tiff to be mosaiked:
            input_image_driver = gdal.Open(equi7_tiff_name, GA_ReadOnly)
            loaded_equi7_tile_data = input_image_driver.GetRasterBand(data_layer_index).ReadAsArray()
            loaded_equi7_tile_quality = input_image_driver.GetRasterBand(quality_layer_index).ReadAsArray()

            # mosaiking (mean where both, use une where just one) it with previous data
            if not "out_mosaiked_data" in locals():
                out_mosaiked_data = loaded_equi7_tile_data
                out_mosaiked_quality = loaded_equi7_tile_quality
            else:
                where_loaded_only_is_not_zero = np.logical_and(out_mosaiked_data == 0, loaded_equi7_tile_data != 0)
                where_both_are_not_zero = np.logical_and(out_mosaiked_data != 0, loaded_equi7_tile_data != 0)

                out_mosaiked_data[where_loaded_only_is_not_zero] = loaded_equi7_tile_data[where_loaded_only_is_not_zero]
                out_mosaiked_data[where_both_are_not_zero] = (
                    loaded_equi7_tile_data[where_both_are_not_zero] + out_mosaiked_data[where_both_are_not_zero]
                ) / 2

                out_mosaiked_quality[where_loaded_only_is_not_zero] = loaded_equi7_tile_quality[
                    where_loaded_only_is_not_zero
                ]
                out_mosaiked_quality[where_both_are_not_zero] = (
                    loaded_equi7_tile_quality[where_both_are_not_zero] + out_mosaiked_quality[where_both_are_not_zero]
                ) / 2

        # write the output final mosaiked tiff of the current equi7 tile:
        out_mosaiked_equi7_name = os.path.join(out_tile_equi7_full_folder, "FH.tif",)

        geotransform = input_image_driver.GetGeoTransform()
        projection = input_image_driver.GetProjection()

        # write the tiff as an equi7
        Nx, Ny = out_mosaiked_data.shape

        driver = gdal.GetDriverByName("GTiff")
        os.makedirs(os.path.dirname(out_mosaiked_equi7_name))
        outdata = driver.Create(out_mosaiked_equi7_name, Ny, Nx, 2, gdal.GDT_Float32)
        outdata.SetGeoTransform(geotransform)
        outdata.SetProjection(projection)
        outdata.GetRasterBand(data_layer_index).WriteArray(out_mosaiked_data)
        outdata.GetRasterBand(quality_layer_index).WriteArray(out_mosaiked_quality)
        outdata.FlushCache()  # saves to disk
        outdata = None

        out_mosaiked_equi7_tiff_names.append(out_mosaiked_equi7_name)

        del out_mosaiked_data

    return out_mosaiked_equi7_tiff_names


def apply_dem_flattening(beta0_calibrated, kz_in, reference_height, master_id, raster_info):

    for pf_name in beta0_calibrated.keys():
        for pol_id in beta0_calibrated[pf_name]:
            beta0_calibrated[pf_name][pol_id] = np.multiply(
                beta0_calibrated[pf_name][pol_id], np.exp(-1j * np.multiply(kz_in[pf_name], reference_height)),
            )

    return beta0_calibrated
