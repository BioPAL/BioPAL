# SPDX-FileCopyrightText: BioPAL <biopal@esa.int>
# SPDX-License-Identifier: MIT

import os
import logging
import numpy as np
import pyproj as proj
from shutil import which
from datetime import datetime
from arepytools.timing.precisedatetime import PreciseDateTime
from arepytools.io.productfolder import ProductFolder
from arepytools.io.metadata import RasterInfo
from equi7grid.equi7grid import Equi7Grid
from biopal.utility.constants import LIGHTSPEED

# The class Task is the template to be inherited for the creation of each BioPAL APP
class Task:
    def __init__(self, configuration_file=None):

        self.configuration_file = configuration_file

    def name(self):
        return self.__class__.__name__

    def _run(self, input_file_xml):
        raise NotImplementedError()

    def run(self, input_file_xml):

        try:

            logging.info("Starting {}".format(self.name()))
            output = self._run(input_file_xml)
            logging.info("Finished {}".format(self.name()))

            return output

        except Exception as e:
            error_msg = "biopal error inside {}: {}".format(self.name(), e)
            logging.error(error_msg, exc_info=True)
            raise RuntimeError(error_msg)


def start_logging(output_folder, L2_product, log_level, app_name=None):
    # CRITICAL 50
    # ERROR 40
    # WARNING 30
    # INFO 20
    # DEBUG 10
    # NOTSET 0

    if log_level == "DEBUG":
        level_to_set = logging.DEBUG
    elif log_level == "INFO":
        level_to_set = logging.INFO
    elif log_level == "WARNING":
        level_to_set = logging.WARNING
    elif log_level == "ERROR":
        level_to_set = logging.ERROR

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    if not app_name:
        log_file_name = os.path.join(output_folder, "biomassL2.log")
    else:
        log_file_name = os.path.join(output_folder, app_name + "_APP.log")

    if os.path.exists(log_file_name):
        raise RuntimeError("output folder {} is not empty.".format(output_folder))

    logging.basicConfig(
        handlers=[logging.FileHandler(log_file_name, mode="w", encoding="utf-8"), logging.StreamHandler(),],
        level=level_to_set,
        format="%(asctime)s - %(levelname)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    logging.getLogger("matplotlib.font_manager").disabled = True

    logging.info(" --BIOMASS L2 Processor-- ")
    if not app_name:
        logging.info("Executing BioPAL full processor, L2 product: {}".format(L2_product))
    else:
        logging.info("Executing {} APP".format(app_name))

    logging.info(" \n")

    return log_file_name


def set_gdal_paths(gdal_path, gdal_environment_path=None):
    if gdal_path is None or not gdal_path:
        gdal_info_path = which("gdalinfo")
        if not gdal_info_path:
            raise RuntimeError("Missing gdalinfo executable")
        gdal_path = os.path.dirname(gdal_info_path)

    if gdal_environment_path is None or not gdal_environment_path:
        gdal_environment_path = os.environ.get("GDAL_DATA")
        if not gdal_environment_path:
            raise RuntimeError("Missing GDAL_DATA environment variable")
    else:
        # Set the enviroment
        os.environ["GDAL_DATA"] = gdal_environment_path

    logging.info("gdal_path: " + gdal_path)
    logging.info("gdal_environment_path: " + gdal_environment_path)

    return gdal_path, gdal_environment_path


def get_data_time_stamp(folder, pf_name):
    # reads a data:
    # it is supposed to contain one or more polarizations (the "SwathInfo" is read to retrive it)
    # it returns a dictionary with keys = polarizations
    # it returns also the dimensions of data

    data_pf_name = os.path.join(folder, pf_name)

    pf = ProductFolder(data_pf_name, "r")

    # prepare the metadata elements
    data_channel_obj = pf.get_channel(0)
    metadata_obj = data_channel_obj.metadata
    metadatachannel_obj = metadata_obj.get_metadata_channels(0)

    # Raster Info
    ri = metadatachannel_obj.get_element("RasterInfo")

    lines_start_utc = str(ri.lines_start)

    return lines_start_utc


def getBinaryNameFromChannelIDX(pf_name, channel_idx):

    if not type(channel_idx) == type(1) or channel_idx < 1:
        error_message = "Input channel_idx should be an integer number greater than zero"
        logging.error(error_message)
        raise ValueError(error_message)

    channel_idx_str = str(channel_idx)
    number_of_zeros = 4 - len(channel_idx_str)
    zeros_str = "0" * number_of_zeros
    binary_name = pf_name + "_" + zeros_str + channel_idx_str

    return binary_name


def decode_unique_acquisition_id_string(unique_acquisition_id_string, output_format="numeric"):

    if output_format == "string":
        global_cycle_idx = unique_acquisition_id_string[3:5]
        heading_deg = unique_acquisition_id_string[8:14]
        rg_swath_idx = unique_acquisition_id_string[20:22]
        rg_sub_swath_idx = unique_acquisition_id_string[30:32]
        az_swath_idx = unique_acquisition_id_string[38:40]
        baseline_idx = unique_acquisition_id_string[45:]

    elif output_format == "numeric":
        global_cycle_idx = int(unique_acquisition_id_string[3:5])
        heading_deg = float(unique_acquisition_id_string[8:14])
        rg_swath_idx = int(unique_acquisition_id_string[20:22])
        rg_sub_swath_idx = int(unique_acquisition_id_string[30:32])
        az_swath_idx = int(unique_acquisition_id_string[38:40])
        baseline_idx = int(unique_acquisition_id_string[45:])

    else:
        raise ValueError(
            'Input output_format = "{}"  not valid, it should be "numeric" or "string"'.format(output_format)
        )

    return (
        global_cycle_idx,
        heading_deg,
        rg_swath_idx,
        rg_sub_swath_idx,
        az_swath_idx,
        baseline_idx,
    )


def check_if_path_exists(path, file_folder_str):
    if not os.path.exists(path):
        log_error_str = file_folder_str + " " + path + " does not exist"

        logging.error(log_error_str)
        raise Exception(log_error_str)


def check_if_geometry_auxiliaries_are_present(
    file_names,
    stack_id,
    acquisitions_pf_names,
    read_ecef=True,
    read_off_nadir=True,
    read_slope=True,
    read_kz=True,
    read_ref_h=True,
    read_dist=False,  # No mandatory, needed only for FH if Spectral Shift Filtering
):

    aux_are_present_flag = True

    if read_ecef:
        pf_name = os.path.basename(file_names.ECEF_grid_file_names[stack_id])
        folder_name = os.path.dirname(file_names.ECEF_grid_file_names[stack_id])
        if not os.path.exists(os.path.join(folder_name, pf_name)):
            aux_are_present_flag = False

    if read_kz:
        pf_name = os.path.basename(file_names.kz_file_names[stack_id])
        folder_name = os.path.dirname(file_names.kz_file_names[stack_id])
        if not os.path.exists(os.path.join(folder_name, pf_name)):
            aux_are_present_flag = False

    if read_off_nadir:
        pf_name = os.path.basename(file_names.off_nadir_angle_file_names[stack_id])
        folder_name = os.path.dirname(file_names.off_nadir_angle_file_names[stack_id])
        if not os.path.exists(os.path.join(folder_name, pf_name)):
            aux_are_present_flag = False

    if read_ref_h:
        pf_name = os.path.basename(file_names.reference_height_file_names[stack_id])
        folder_name = os.path.dirname(file_names.reference_height_file_names[stack_id])
        if not os.path.exists(os.path.join(folder_name, pf_name)):
            aux_are_present_flag = False

    if read_dist:
        pf_name = os.path.basename(file_names.slant_range_distances_file_names[stack_id])
        folder_name = os.path.dirname(file_names.slant_range_distances_file_names[stack_id])
        if not os.path.exists(os.path.join(folder_name, pf_name)):
            warning_message = " Spectral Shift Filtering is enabled, the Distances Auxiliary product is mandatory in this case but it is missing."
            logging.warning(warning_message)
            aux_are_present_flag = False

    if read_slope:
        pf_name = os.path.basename(file_names.slope_file_names[stack_id])
        folder_name = os.path.dirname(file_names.slope_file_names[stack_id])
        if not os.path.exists(os.path.join(folder_name, pf_name)):
            aux_are_present_flag = False

    return aux_are_present_flag


def choose_equi7_sampling(product_resolution, geographic_grid_sampling):

    if geographic_grid_sampling <= product_resolution / 2:
        equi7_sampling = geographic_grid_sampling
    else:
        equi7_sampling = product_resolution / 2
        warning_message = "user rquested Geographic Grid Sampling is {} [m]; it cannot be greater than (product_resolution)/2 = {} [m] \n".format(
            geographic_grid_sampling, product_resolution / 2
        )
        logging.warning(warning_message)

    _static_sampling_sorted = np.sort(Equi7Grid._static_sampling)
    indexes_less_or_equal = np.where(np.less_equal(_static_sampling_sorted, equi7_sampling))

    if len(indexes_less_or_equal[0]):
        equi7_sampling_quantized = _static_sampling_sorted[indexes_less_or_equal[0][-1]]
    else:
        equi7_sampling_quantized = _static_sampling_sorted[0]

    return equi7_sampling_quantized


def evaluate_estimation_quality_matrix(data_in_shape):
    # Placemark for the quality estimation to be defined

    logging.warning(
        "The quality estimation of the product is still to be defined: zeros will be placed in the quality layer "
    )
    quality_matrix = np.zeros(data_in_shape)

    return quality_matrix


def epsg_in_to_epsg_out(xx, yy, zz, epsg_code_in, epsg_code_out):
    p_in = proj.Proj("+init={}".format(epsg_code_in))
    p_out = proj.Proj("+init={}".format(epsg_code_out))
    xx_transf, yy_transf, zz_transf = proj.transform(p_in, p_out, xx, yy, zz)

    return xx_transf, yy_transf, zz_transf


def convert_rasterinfo_meters_to_seconds(ri_meters, sensor_velocity):

    ri_seconds = RasterInfo(
        ri_meters.lines,
        ri_meters.samples,
        ri_meters.cell_type,
        ri_meters.file_name,
        ri_meters.header_offset_bytes,
        ri_meters.row_prefix_bytes,
        ri_meters.byte_order,
    )

    lines_step_s = ri_meters.lines_step / sensor_velocity
    lines_step_unit = "s"
    samples_start_s = ri_meters.samples_start / LIGHTSPEED * 2
    samples_start_unit = "s"
    samples_step_s = ri_meters.samples_step / LIGHTSPEED * 2
    samples_step_unit = "s"

    ri_seconds.set_lines_axis(ri_meters.lines_start, ri_meters.lines_start_unit, lines_step_s, lines_step_unit)
    ri_seconds.set_samples_axis(samples_start_s, samples_start_unit, samples_step_s, samples_step_unit)

    return ri_seconds


def convert_rasterinfo_seconds_to_meters(ri_seconds, sensor_velocity):

    ri_meters = RasterInfo(
        ri_seconds.lines,
        ri_seconds.samples,
        ri_seconds.cell_type,
        ri_seconds.file_name,
        ri_seconds.header_offset_bytes,
        ri_seconds.row_prefix_bytes,
        ri_seconds.byte_order,
    )

    lines_step_m = ri_meters.lines_step * sensor_velocity
    lines_step_unit = "m"
    samples_start_m = ri_meters.samples_start * LIGHTSPEED / 2
    samples_start_unit = "m"
    samples_step_m = ri_meters.samples_step * LIGHTSPEED / 2
    samples_step_unit = "m"

    ri_meters.set_lines_axis(
        ri_seconds.lines_start, ri_seconds.lines_start_unit, lines_step_m, lines_step_unit,
    )
    ri_meters.set_samples_axis(samples_start_m, samples_start_unit, samples_step_m, samples_step_unit)

    return ri_meters


def save_breakpoints(output_folder, breakpoint_names_list, breakpoint_data_list):

    if len(breakpoint_names_list) != len(breakpoint_data_list):
        raise IndexError("Inputs should have same shape")

    for idx, file_name in enumerate(breakpoint_names_list):
        np.save(os.path.join(output_folder, file_name), breakpoint_data_list[idx])


def format_folder_name():
    # date string format is: BIOMASS_L2_YYYYMMDDTHHMMSS
    now = datetime.now()
    year_str = str(now.year)
    month_str = str(now.month).rjust(2, "0")
    day_str = str(now.day).rjust(2, "0")
    h_str = str(now.hour).rjust(2, "0")
    m_str = str(now.minute).rjust(2, "0")
    s_str = str(now.second).rjust(2, "0")
    current_date_string = "BIOMASS_L2_" + year_str + month_str + day_str + "T" + h_str + m_str + s_str

    return current_date_string


def check_fnf_folder_format(fnf_mask_catalogue_folder):

    fnf_content_folders = os.listdir(fnf_mask_catalogue_folder)
    if not fnf_content_folders:
        error_msg = "input fnf mask is mandatory"
        logging.error(error_msg)
        raise Exception(error_msg)

    found_tandemx = False
    found_equi7 = False
    for fnf_content_folder in fnf_content_folders:
        if "TDM_FNF_" in fnf_content_folder:
            found_tandemx = True

        elif "EQUI7_" in fnf_content_folder:
            found_equi7 = True

        else:
            error_msg = 'Cannot recognize input FNF mask format: valid formats are EQUI7 (an "EQUI7_*/" folder containing equi7 tiles sub-folders, example: "EQUI7_AF025M/E042N048T3) or TANDEM-X (collection of "TDM_FNF_*/" folders '
            logging.error(error_msg)
            raise Exception(error_msg)

    if found_tandemx and found_equi7:
        error_msg = " Found bowth TANDEM-X and EQUI7 format FNF masks, this is not supported."
        logging.error(error_msg)
        raise Exception(error_msg)
    elif found_tandemx:
        return "TANDEM-X"
    elif found_equi7:
        return "EQUI7"


def check_cal_format(reference_agb_folder):

    cal_format = "RASTER"

    return cal_format


def get_min_time_stamp_repository(L1c_repository, stack_composition):

    for idx, (unique_stack_id, unique_acq_pf_names) in enumerate(stack_composition.items()):

        if idx == 0:
            time_tag_mjd_initial = get_min_time_stamp_repository_core(L1c_repository, unique_acq_pf_names)
        else:
            time_tag_mjd_initial = min(
                time_tag_mjd_initial, get_min_time_stamp_repository_core(L1c_repository, unique_acq_pf_names),
            )

    time_tag_mjd_initial = PreciseDateTime().set_from_utc_string(time_tag_mjd_initial)

    return time_tag_mjd_initial


def get_min_time_stamp_repository_core(L1c_repository, acquisitions_pf_names):
    for idx, pf_name in enumerate(acquisitions_pf_names):

        if idx == 0:
            min_time_stamp_mjd = get_data_time_stamp(L1c_repository, pf_name)
        else:
            min_time_stamp_mjd = min(min_time_stamp_mjd, get_data_time_stamp(L1c_repository, pf_name))

    return min_time_stamp_mjd


def check_equi7_mask_coverage(data_equi7_fnames, equi7_fnf_mask_fnames):
    equi7_tile_names = []
    for fname_curr in data_equi7_fnames:
        equi7_tile_names.append(os.path.basename(os.path.dirname(fname_curr)))
    input_mask_equi7_tile_names = []
    for m_name_curr in equi7_fnf_mask_fnames:
        input_mask_equi7_tile_names.append(os.path.basename(os.path.dirname(m_name_curr)))
    for equi7_tile_name in equi7_tile_names:
        if not equi7_tile_name in input_mask_equi7_tile_names:
            error_msg = 'The auxiliary "ForestMask" does not cover the acquisitions, missing Equi7 tile "{}"'.format(
                equi7_tile_name
            )
            logging.error(error_msg)
            raise ValueError(error_msg)


def resolution_heading_correction(resolution, heading_deg):

    if heading_deg < 0 or heading_deg > 360:
        error_str = "Stack heading not valid: {}  [deg]".format(heading_deg)
        logging.error(error_str)
        raise ValueError(error_str)

    if heading_deg >= 0 and heading_deg <= 45:
        correction_factor = np.cos(np.deg2rad(heading_deg))

    elif heading_deg > 45 and heading_deg <= 135:
        correction_factor = np.sin(np.deg2rad(heading_deg))

    elif heading_deg > 135 and heading_deg <= 225:
        correction_factor = np.abs(np.cos(np.deg2rad(heading_deg)))

    elif heading_deg > 225 and heading_deg <= 315:
        correction_factor = np.abs(np.sin(np.deg2rad(heading_deg)))

    elif heading_deg > 315:
        correction_factor = np.cos(np.deg2rad(heading_deg))

    resolution = resolution / correction_factor

    return resolution


def get_equi7_fnf_tiff_names(fnf_mask_catalogue_folder):

    equi7_tiles_subfolders = []
    e7g_type = ""
    for idx, equi7_subgrid_name in enumerate(os.listdir(fnf_mask_catalogue_folder)):

        if idx == 1:
            logging.warning(
                "Found more than one sub-grid in input Equi7 format FNF mask: prototype supports just one grid at a time, only {} sub-grid will be read".format(
                    equi7_subgrid_name
                )
            )
            break

        equi7_folder = os.path.join(fnf_mask_catalogue_folder, equi7_subgrid_name)
        temp_list = os.listdir(equi7_folder)
        if e7g_type == "":
            e7g_type = temp_list[0][-2:]
        if e7g_type in temp_list[0][-2:]:
            for folder_curr in temp_list:
                equi7_tiles_subfolders.append(os.path.join(equi7_folder, folder_curr))
        elif not e7g_type == temp_list[0][-2:]:
            logging.warning(
                'Auxiliary forest mask multiple tiles sampling are not supported, only "'
                + e7g_type
                + '" tile sampling have been read, "'
                + temp_list[0][-2:]
                + '" skipped instead'
            )

    equi7_fnf_mask_tiff_names_list = []
    for equi7_tile_folder in equi7_tiles_subfolders:

        tiff_fname = os.listdir(equi7_tile_folder)[0]

        equi7_fnf_mask_tiff_names_list.append(os.path.join(equi7_tile_folder, tiff_fname))

    return equi7_fnf_mask_tiff_names_list, equi7_subgrid_name


def get_raster_cal_names(reference_agb_folder):

    tiff_names = os.listdir(reference_agb_folder)

    cal_names = []
    for tiff_name in tiff_names:
        if ".tif" in tiff_name or ".tiff" in tiff_name:
            cal_names.append(os.path.join(reference_agb_folder, tiff_name))

    return cal_names


def get_foss_cal_names(reference_agb_folder):

    tiff_names = os.listdir(reference_agb_folder)

    cal_names = []
    for tiff_name in tiff_names:
        if ".tif" in tiff_name or ".tiff" in tiff_name:
            cal_names.append(os.path.join(reference_agb_folder, tiff_name))

    return cal_names


def collect_stacks_to_be_merged(stack_composition):
    """
    The stack_composition contains the IDs of all the stacks and
    for each of the stacks contains all the acquisition IDs.

    This function collects togeter all the stacks which have all in common except the heading:
    the resulting stacks_to_merge_dict is a dictionary containing
        as keys,   all the unique_merged_stack_id (stack IDs with missing headings)
        as values, the list of the stack_ids to be merged (same stacks, all headings)

    """
    stacks_to_merge_dict = {}
    for unique_stack_id in stack_composition.keys():

        # 1) read the acquisition ID of first acquisition in the current stack
        (
            global_cycle_idx,
            heading_deg,
            rg_swath_idx,
            rg_sub_swath_idx,
            az_swath_idx,
            baseline_idx,
        ) = decode_unique_acquisition_id_string(stack_composition[unique_stack_id][0], output_format="string")

        # 2) create a new ID which is a stack_id without the heading
        unique_merged_stack_id = (
            "GC_" + global_cycle_idx + "_RGSW_" + rg_swath_idx + "_RGSBSW_" + rg_sub_swath_idx + "_AZSW_" + az_swath_idx
        )

        if unique_merged_stack_id in stacks_to_merge_dict.keys():
            stacks_to_merge_dict[unique_merged_stack_id].append(unique_stack_id)
        else:
            stacks_to_merge_dict[unique_merged_stack_id] = [unique_stack_id]

    return stacks_to_merge_dict


def radiometric_correction_beta_to_sigma(beta_nought, theta):

    return beta_nought * np.sin(theta)

