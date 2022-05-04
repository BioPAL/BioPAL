# SPDX-FileCopyrightText: BioPAL <biopal@esa.int>
# SPDX-License-Identifier: MIT

import os
import ast
import logging
import numpy as np
import xml.etree.ElementTree as ET
from xml.etree.ElementTree import ElementTree, Element, SubElement
from collections import namedtuple
from namedlist import namedlist
from arepytools.timing.precisedatetime import PreciseDateTime


###############################################################################
# structures (namedtuples) for read and write inputs and configuration files ##
_XML_VERSION = "0.1"

geographic_boundaries_latlon = namedlist(
    "geographic_boundaries",
    "lon_min \
     lon_max \
     lat_min \
     lat_max",
)
geographic_boundaries_northeast = namedlist(
    "geographic_boundaries",
    "east_min \
     east_max \
     north_min \
     north_max",
)

# initialize output structure fields
raster_info = namedtuple(
    "raster_info",
    "num_samples \
 num_lines \
 pixel_spacing_slant_rg \
 pixel_spacing_az \
 resolution_m_slant_rg \
 resolution_m_az \
 carrier_frequency_hz \
 range_bandwidth_hz \
 lines_start_utc",
)

# Main Input parameters structure (namedtuple)
input_params = namedlist(
    "input_params",
    "L2_product \
 output_specification \
 dataset_query \
 stack_based_processing \
 core_processing_agb \
 core_processing_fh \
 core_processing_fd \
 core_processing_tomo_fh \
 core_processing_fnf",
)

output_specification = namedlist(
    "output_specification",
    "output_folder \
    geographic_grid_sampling",
)
dataset_query = namedtuple(
    "dataset_query",
    "L1C_repository \
    L1C_aux_repository \
    L1C_date \
    boundaries \
    stacks_to_find",
)
stack_based_processing = namedtuple(
    "stack_based_processing",
    "stack_composition \
 ECEF_grid_file_names \
 kz_file_names \
 slant_range_distances_file_names \
 off_nadir_angle_file_names \
 slope_file_names \
 reference_height_file_names \
 calibration_screens_file_names \
 dem_folder \
 reference_agb_folder \
 forest_mask_catalogue_folder \
 system_decorrelation_fun_folder \
 average_covariance_folder \
 forest_height_folder \
 geographic_boundaries \
 geographic_boundaries_per_stack",
)
core_processing_agb = namedtuple(
    "core_processing_agb",
    "lut_cal \
    lut_fnf \
    lut_stacks",
)
lut = namedtuple(
    "lut",
    "paths \
     boundaries \
     progressive",
)
core_processing_fh = namedtuple("core_proc_fh", "data_equi7_fnames mask_equi7_fnames")
core_processing_tomo_fh = namedtuple("core_proc_tomo_fh", "data_equi7_fnames mask_equi7_fnames")
core_processing_fd = namedtuple("core_proc_fd", "cycles_composition",)
core_processing_fnf = namedtuple("core_proc_fnf", "pforest_equi7_fnames mask_equi7_fnames")
# main_input_params "L1c_date" sub-fields:
L1c_date = namedtuple(
    "L1c_date_params",
    "start \
     stop",
)
# main_input_params "L1c_date" sub-fields:
geographic_boundary = namedtuple("geographic_boundary_params", "point")

# Configuration parameters structure (namedtuple)
conf_params = namedtuple(
    "conf_params",
    "gdal \
     processing_flags \
     ground_cancellation \
     estimate_agb \
     estimate_fh \
     estimate_tomo_fh \
     change_detection_fd \
     estimate_fnf",
)
conf_gdal = namedtuple(
    "conf_gdal",
    "gdal_path \
     gdal_environment_path",
)
conf_flags = namedtuple(
    "conf_flags",
    "enable_resampling \
     compute_geometry \
	 apply_calibration_screen \
	 DEM_flattening \
	 multilook_heading_correction \
     input_data_type \
	 save_breakpoints \
	 delete_temporary_files",
)

conf_ground_canc = namedtuple(
    "conf_ground_canc",
    "multi_master_flag \
     enhanced_forest_height \
     equalization_flag",
)
# configuration_params "agb" sub-fields:
conf_agb_est = namedtuple(
    "conf_agb_est",
    "product_resolution \
    forest_class_observable_name \
	transfer_function_name \
	number_of_tests \
	fraction_of_roi_per_test \
	fraction_of_cal_per_test \
	add_variability_on_cal_data \
	intermediate_ground_averaging \
	distance_sampling_area \
	parameter_block_size \
	distance_parameter_block \
	min_number_of_rois \
	min_number_of_rois_per_stack \
	min_number_of_cals_per_test \
	min_number_of_rois_per_test \
    estimation_valid_values_limits \
    residual_function",
)
# agb "residual_function" sub-fields:
conf_residual_function = namedtuple(
    "conf_residual_function",
    "formula_terms \
    formula_parameters \
    formula_observables",
)
# agb "formula_terms" sub-fields:
formula_terms = namedtuple(
    "formula_terms",
    "name \
     string \
     formula_weights",
)
formula_weights = namedtuple(
    "formula_weights",
    "step1 \
     step2 \
     step3",
)
# agb "formula_parameters" sub-fields:
formula_parameters = namedtuple(
    "formula_parameters",
    "name \
     save_as_map \
     transform \
     limits \
     limit_units \
     associated_observable_name \
     parameter_variabilities",
)
parameter_variabilities = namedtuple(
    "parameter_variabilities",
    "samples \
     forest_classes \
     stacks \
     global_cycles \
     headings \
     swaths \
     subswath \
     azimuth_images",
)
# agb "formula_observables" sub-fields:
formula_observables = namedtuple(
    "formula_observables",
    "name \
     is_required \
     source_paths \
     source_unit \
     source_resolution \
     limit_units \
     limits \
     transform \
     averaging_method",
)
# min_max is used for parameters limits and observables ranges
min_max = namedlist(
    "min_max",
    "min \
     max",
)

# model_parameters "triplet_params" sub-fields:
triplet_params = namedtuple("triplet_params", "agb_model_scaling agb_model_exponent agb_cosine_exponent")

# configuration_params "FH" sub-fields:
conf_fh_est = namedtuple(
    "conf_fh_est",
    "spectral_shift_filtering \
     product_resolution \
     kz_thresholds \
     model_parameters \
     median_factor",
)
# FH "model_parameters" sub-fields:
conf_fh_model_params = namedtuple(
    "conf_fh_model_params",
    "estimation_valid_values_limits \
     maximum_height \
     number_of_extinction_value \
     number_of_ground_volume_ratio_value \
     number_of_temporal_decorrelation_value",
)

# configuration_params "FD" sub-fields:
conf_fd_est = namedtuple(
    "conf_fd_est",
    "product_resolution \
    confidence_level",
)
# configuration_params "TOMO FH" sub-fields:
conf_tomo_fh_est = namedtuple(
    "conf_tomo_fh_est",
    "product_resolution \
     vertical_range \
     estimation_valid_values_limits \
     enable_super_resolution \
     regularization_noise_factor \
     power_threshold \
     median_factor",
)

# TOMO and TOMO_FH "vertical range" sub-fields:
vertical_range_params = namedtuple(
    "vertical_range_params",
    "maximum_height \
     minimum_height \
     sampling",
)
# configuration_params "FNF" sub-fields:
conf_fnf_est = namedtuple(
    "conf_fnf_est",
    "product_resolution \
    logreg_coeffs_fname",
)
###############################################################################


# BioPAL Input File parser and writer
def write_input_file(input_params_obj, input_file_xml):
    """
    Write the input file:
        
    the input file is divided in many sections;
        - L2_product (contains the flag to enable specific chains, AGB, FH, FD or TOMO_FH)
        - output_specification
        - dataset_query
        - stack_based_processing (internally the presence or not of specific 
                                  elements vary on the chain, AGB,FH,FD,TOMO_FH)
        - core_processing_agb
        - core_processing_fh
        - core_processing_fd
        - core_processing_tomo_fh
        - core_processing_fnf
    
    Only the sections avalable into the input_params_obj structure will be written 
    to the file: not all the APPs requires all the sections. 
    
    INPUT:  input_param_obj, a structure, as generated by read_inpuit_file()
            path of the xml input file to be written
        
    See Also:
        read_input_file(), which takes in input the xml path and returns input_param_obj
    """

    root_item = Element("BioPALInput")
    root_item.set("version", _XML_VERSION)

    # write L2_product section
    if input_params_obj.L2_product:
        write_L2_product_section(root_item, input_params_obj.L2_product)

    # write output_specification section
    if input_params_obj.output_specification:
        write_output_specification_section(root_item, input_params_obj.output_specification)

    # write dataset_query section
    if input_params_obj.dataset_query:
        write_dataset_query_section(root_item, input_params_obj.dataset_query)

    # write stack_based_processing section
    if input_params_obj.stack_based_processing:
        write_stack_based_processing_section(root_item, input_params_obj.stack_based_processing)

    # write stack_based_processing section
    if input_params_obj.core_processing_agb:
        write_core_processing_agb_section(root_item, input_params_obj.core_processing_agb)

    # write stack_based_processing section
    if input_params_obj.core_processing_fh:
        write_core_processing_fh_section(root_item, input_params_obj.core_processing_fh)

    # write stack_based_processing section
    if input_params_obj.core_processing_fd:
        write_core_processing_fd_section(root_item, input_params_obj.core_processing_fd)

    # write stack_based_processing section
    if input_params_obj.core_processing_tomo_fh:
        write_core_processing_tomo_fh_section(root_item, input_params_obj.core_processing_tomo_fh)
        
    # write stack_based_processing section
    if input_params_obj.core_processing_fnf:
        write_core_processing_fnf_section(root_item, input_params_obj.core_processing_fnf)

    # write to file
    output_dir = os.path.dirname(input_file_xml)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    ElementTree_indent(root_item)
    tree = ElementTree(root_item)
    tree.write(open(input_file_xml, "w"), encoding="unicode")


def write_L2_product_section(father_item, L2_product_string):

    L2_product_item = SubElement(father_item, "L2_product")
    L2_product_item.set("options", "AGB, FH, FD, TOMO_FH")
    L2_product_item.text = L2_product_string


def write_output_specification_section(father_item, output_spec_obj):

    out_spec_item = SubElement(father_item, "output_specification")

    out_folder_item = SubElement(out_spec_item, "output_folder")
    out_folder_item.text = output_spec_obj.output_folder
    geo_samp_item = SubElement(out_spec_item, "geographic_grid_sampling")
    geo_samp_item.text = str(output_spec_obj.geographic_grid_sampling)
    geo_samp_item.set("unit", "m")


def write_dataset_query_section(father_item, dataset_query_obj):

    dataset_query_item = SubElement(father_item, "dataset_query")

    L1C_repository_item = SubElement(dataset_query_item, "L1C_repository")
    L1C_repository_item.text = dataset_query_obj.L1C_repository

    L1C_date_start_item = SubElement(dataset_query_item, "L1C_date")
    L1C_date_start_item.set("value", "start")
    L1C_date_start_item.set("unit", "Utc")
    L1C_date_start_item.text = str(dataset_query_obj.L1C_date.start)

    L1C_date_stop_item = SubElement(dataset_query_item, "L1C_date")
    L1C_date_stop_item.set("value", "stop")
    L1C_date_stop_item.set("unit", "Utc")
    L1C_date_stop_item.text = str(dataset_query_obj.L1C_date.stop)

    aux_folder_item = SubElement(dataset_query_item, "auxiliary_products_folder")
    aux_folder_item.text = dataset_query_obj.L1C_aux_repository

    geo_bound_poly_item = SubElement(dataset_query_item, "geographic_boundaries_polygon")
    geo_bound_poly_item.set("unit", "deg")
    for point in dataset_query_obj.boundaries.point:
        point_item = SubElement(geo_bound_poly_item, "point")
        lat_item = SubElement(point_item, "latitude")
        lat_item.text = str(point["Latitude"])
        lon_item = SubElement(point_item, "longitude")
        lon_item.text = str(point["Longitude"])

    if dataset_query_obj.stacks_to_find:
        stacks_to_find_item = SubElement(dataset_query_item, "stacks_to_find")
        stacks_to_find_item.text = dataset_query_obj.stacks_to_find


def write_stack_based_processing_section(father_item, stack_based_proc_obj):

    sb = stack_based_proc_obj

    stack_based_item = SubElement(father_item, "stack_based_processing")

    write_boundaries_latlon_section(stack_based_item, sb.geographic_boundaries)

    L1C_product_list_item = SubElement(stack_based_item, "L1C_product_list")

    for stack_id, stack_composition_list in sb.stack_composition.items():

        L1C_product_item = SubElement(L1C_product_list_item, "L1C_product")
        L1C_product_item.set("unique_stack_id", stack_id)

        for acq_name in stack_composition_list:
            acquisition_item = SubElement(L1C_product_item, "acquisition")
            acquisition_item.text = acq_name

        write_boundaries_latlon_section(L1C_product_item, sb.geographic_boundaries_per_stack[stack_id])

        auxiliaries_item = SubElement(L1C_product_item, "auxiliaries")
        geometry_item = SubElement(auxiliaries_item, "geometry")

        ECEF_item = SubElement(geometry_item, "ECEF_grid")
        ECEF_item.text = sb.ECEF_grid_file_names[stack_id]
        vw_item = SubElement(geometry_item, "vertical_wavenumber")
        vw_item.text = sb.kz_file_names[stack_id]
        rd_item = SubElement(geometry_item, "radar_distances")
        rd_item.text = sb.slant_range_distances_file_names[stack_id]
        ona_item = SubElement(geometry_item, "off_nadir_angles")
        ona_item.text = sb.off_nadir_angle_file_names[stack_id]
        slope_item = SubElement(geometry_item, "slope")
        slope_item.text = sb.slope_file_names[stack_id]
        rh_item = SubElement(geometry_item, "reference_height")
        rh_item.text = sb.reference_height_file_names[stack_id]

        if sb.calibration_screens_file_names:
            cal_screens_item = SubElement(auxiliaries_item, "calibration_screens")
            cal_screens_item.text = sb.calibration_screens_file_names[stack_id]

    dem_item = SubElement(stack_based_item, "DEM")
    dem_item.text = sb.dem_folder

    if sb.reference_agb_folder:
        ref_agb_item = SubElement(stack_based_item, "reference_agb")
        ref_agb_item.text = sb.reference_agb_folder
    if sb.forest_height_folder:
        fh_folder_item = SubElement(stack_based_item, "forest_height")
        fh_folder_item.text = sb.forest_height_folder
    if sb.forest_mask_catalogue_folder:
        forest_mask_item = SubElement(stack_based_item, "forest_mask")
        forest_mask_item.text = sb.forest_mask_catalogue_folder
    if sb.system_decorrelation_fun_folder:
        system_dec_fun_item = SubElement(stack_based_item, "system_decorrelation_function")
        system_dec_fun_item.text = sb.system_decorrelation_fun_folder


def write_core_processing_agb_section(father_item, core_proc_agb_obj):

    core_agb_item = SubElement(father_item, "core_processing_agb")

    write_lut_section(core_agb_item, core_proc_agb_obj.lut_cal, "cal")
    write_lut_section(core_agb_item, core_proc_agb_obj.lut_fnf, "fnf")
    write_lut_section(core_agb_item, core_proc_agb_obj.lut_stacks, "stacks")


def write_core_processing_fh_section(father_item, core_proc_fh_obj, tomo_fh_flag=False):

    if tomo_fh_flag:
        core_fh_item = SubElement(father_item, "core_processing_tomo_fh")
    else:
        core_fh_item = SubElement(father_item, "core_processing_fh")

    equi7_products_paths_item = SubElement(core_fh_item, "equi7_products_paths")

    for stack_id, path_data_fnames in core_proc_fh_obj.data_equi7_fnames.items():

        equi7_product_item = SubElement(equi7_products_paths_item, "equi7_product")
        equi7_product_item.set("unique_stack_id", stack_id)

        # each stack contains n-equi7 tiles
        for tile_idx, path_mask_tile in enumerate(core_proc_fh_obj.mask_equi7_fnames[stack_id]):

            path_data_tile = path_data_fnames[tile_idx]

            tile_item = SubElement(equi7_product_item, "equi7_tile")

            path_data_item = SubElement(tile_item, "path_data")
            path_data_item.text = path_data_tile

            path_mask_item = SubElement(tile_item, "path_mask")
            path_mask_item.text = path_mask_tile

def write_core_processing_fnf_section(father_item, core_proc_fnf_obj):
    
    core_fnf_item = SubElement(father_item, "core_processing_fnf")

    equi7_products_paths_item = SubElement(core_fnf_item, "equi7_products_paths")

    for stack_id, path_data_fnames in core_proc_fnf_obj.pforest_equi7_fnames.items():
        
        equi7_product_item = SubElement(equi7_products_paths_item, "equi7_product")
        equi7_product_item.set("unique_stack_id", stack_id)
        
        # each stack contains n-equi7 tiles
        for tile_idx, path_mask_tile in enumerate(core_proc_fnf_obj.mask_equi7_fnames[stack_id]):
            
            path_data_tile = path_data_fnames[tile_idx]
            
            tile_item = SubElement(equi7_product_item, "equi7_tile")
            
            path_data_item = SubElement(tile_item, "path_data")
            path_data_item.text = path_data_tile
            
            path_mask_item = SubElement(tile_item, "path_mask")
            path_mask_item.text = path_mask_tile


def write_core_processing_fd_section(father_item, core_proc_fd_obj):

    core_fd_item = SubElement(father_item, "core_processing_fd")

    cycles_composition_item = SubElement(core_fd_item, "cycles_composition")
    # dict which groupes together all the differente global cycles for stacks with same nominal geometry
    # cycles_composition[nominal_geometry_stack_id][global_cycle_number]

    for stack_idx, (nominal_geometry_stack_id, global_cycle_dict) in enumerate(
        core_proc_fd_obj.cycles_composition.items()
    ):

        nominal_geometry_item = SubElement(cycles_composition_item, "nominal_geometry")
        nominal_geometry_item.set("stack_id", nominal_geometry_stack_id)

        for cycle_number, acquisitions_list in global_cycle_dict.items():

            global_cycle_item = SubElement(nominal_geometry_item, "global_cycle")
            global_cycle_item.set("cycle_number", str(cycle_number))

            for acquisition_id in acquisitions_list:
                acquisition_id_item = SubElement(global_cycle_item, "acquisition")
                acquisition_id_item.text = acquisition_id


def write_core_processing_tomo_fh_section(father_item, core_proc_fh_obj):

    # core processing section for FH and TOMO FH are the same
    write_core_processing_fh_section(father_item, core_proc_fh_obj, tomo_fh_flag=True)


def parse_input_file(input_file_xml):
    """
    Parse the input file:
        
    the input file is divided in many sections;
    the three mandatory sections to launch the main BioPAL APP or the dataset_query APP are:
        - L2_product (contains the flag to enable specific chains, AGB, FH, FD or TOMO_FH)
        - output_specification
        - dataset_query
    Other sections available are:
        - stack_based_processing (internally the presence of specific 
                                  elements vary on the chain, AGB,FH,FD,TOMO_FH)
        - core_processing_agb
        - core_processing_fh (place holder, to be filled in future version if needed)
        - core_processing_fd (place holder, to be filled in future version if needed)
        - core_processing_tomo_fh (place holder, to be filled in future version if needed)
    
    Some APPs can update the input file, adding new sections:
        i.e. the stack based processing APP of the AGB will write the new 
        section "core_processing_agb"  to be passed to the AGB core APP 
        
    In general each section of the input file is considered optional in this parser:
    each APP will check if the needed sections are present and compliant or not.
    
    INPUT: path of the xml input file
    OUTPUT: input_param_obj, a structure containing all the fields
    
    See Also:
        write_input_file(), which takes in input the input_param_obj and writes to xml
    """

    tree = ET.parse(input_file_xml)
    root = tree.getroot()

    # L2_product section
    L2_product_string = parse_L2_product_section(root)

    # output_specification section
    output_spec_obj = parse_output_specification_section(root)

    # dataset_query section
    dataset_query_obj = parse_dataset_query_section(root)

    # stack_based_processing section
    stack_based_proc_obj = parse_stack_based_proc_section(root)

    # core_processing_agb section
    core_proc_agb_obj = parse_core_proc_agb_section(root)

    # core_processing_fh section
    core_proc_fh_obj = parse_core_proc_fh_section(root)

    # core_processing_fd section
    core_proc_fd_obj = parse_core_proc_fd_section(root)

    # core_processing_tomo_fh section
    core_proc_tomo_fh_obj = parse_core_proc_tomo_fh_section(root)
    
    # core_processing_fnf section
    core_proc_fnf_obj = parse_core_proc_fnf_section(root)

    # output structures filling
    input_params_obj = input_params(
        L2_product_string,
        output_spec_obj,
        dataset_query_obj,
        stack_based_proc_obj,
        core_proc_agb_obj,
        core_proc_fh_obj,
        core_proc_fd_obj,
        core_proc_tomo_fh_obj,
        core_proc_fnf_obj
    )

    return input_params_obj


def parse_L2_product_section(root):

    L2Product_Item = root.find("L2_product")
    if not L2Product_Item is None:
        L2_product_string = L2Product_Item.text
        if not (
            L2_product_string == "AGB"
            or L2_product_string == "FH"
            or L2_product_string == "FD"
            or L2_product_string == "TOMO_FH"
            or L2_product_string == "FNF"
        ):

            error_message = [L2_product_string + " is an invalid L2_product, possible values are AGB, FH, FD, TOMO_FH, FNF"]
            logging.error(error_message)
            raise ValueError(error_message)
    else:
        L2_product_string = None

    return L2_product_string


def parse_output_specification_section(root):

    output_specification_item = root.find("output_specification")

    if output_specification_item:
        output_folder = output_specification_item.find("output_folder").text
        geographic_grid_sampling = float(output_specification_item.find("geographic_grid_sampling").text)

        output_spec_obj = output_specification(output_folder, geographic_grid_sampling)
    else:
        output_spec_obj = None

    return output_spec_obj


def parse_dataset_query_section(root):

    dataset_query_item = root.find("dataset_query")

    if dataset_query_item:

        L1c_repository = dataset_query_item.find("L1C_repository").text
        if not os.path.exists(L1c_repository):
            error_message = "Main input file: the specified L1C_repository folder does not exist "
            logging.error(error_message)
            raise ValueError(error_message)

        L1cDates = dataset_query_item.findall("L1C_date")
        if len(L1cDates) != 2:
            error_message = (
                'Main input file: there should be #2 L1C dates (one with value="start" and one with value="stop"'
            )
            logging.error(error_message)
            raise ValueError(error_message)
        for date_UTC_Item in L1cDates:

            L1cDate = date_UTC_Item.text
            attrib = date_UTC_Item.attrib["value"]

            if "start" in attrib:
                L1cDate_start = PreciseDateTime().set_from_utc_string(L1cDate)
            elif "stop" in attrib:
                L1cDate_stop = PreciseDateTime().set_from_utc_string(L1cDate)
            else:
                logging.error(error_message)
                raise ValueError(
                    'Main input file: L1cDate "value" attribute should be "start" or "stop": "'
                    + attrib
                    + '" is not supported'
                )

        if L1cDate_stop < L1cDate_start:
            error_message = 'Main input file: L1cDate "start" value cannot be greater of "stop" value'
            logging.error(error_message)
            raise ValueError(error_message)

        L1c_date_struct = L1c_date(L1cDate_start, L1cDate_stop)

        L1c_aux_data_repository = dataset_query_item.find("auxiliary_products_folder").text
        if not os.path.exists(L1c_aux_data_repository):
            error_message = "Main input file: the specified auxiliary_products_folder does not exist: "
            logging.error(error_message)
            raise ValueError(error_message)

        dem_folder_found = False
        aux_list = os.listdir(L1c_aux_data_repository)
        for aux_name in aux_list:
            if aux_name == "DEM":
                dem_folder_found = True
                break
        if not dem_folder_found:
            error_message = (
                'Main input file: the specified auxiliary_products_folder should contain AT LEAST the sub-folder "DEM"'
            )
            logging.error(error_message)
            raise ValueError(error_message)

        geographic_boundary_struct = parse_polygon_boundaries_section(dataset_query_item)

        # the following is an optional field, for debug purposes only
        stacks_to_find_list = []
        Stacks_to_find_Item = dataset_query_item.find("stacks_to_find")
        if Stacks_to_find_Item != None:
            stacks_to_find = Stacks_to_find_Item.findall("stack_id")
            if len(stacks_to_find):
                for stack_curr in stacks_to_find:
                    stacks_to_find_list.append(stack_curr.text)
                    if stack_curr.text == None:
                        error_message = 'Main input file: "Stacks to find" optional input cannot have empty "Stack_id" fields: fill them or totally erase the "Stacks to find" element from input file'
                        logging.error(error_message)
                        raise ValueError(error_message)

        dataset_query_obj = dataset_query(
            L1c_repository, L1c_aux_data_repository, L1c_date_struct, geographic_boundary_struct, stacks_to_find_list
        )

    else:

        dataset_query_obj = None

    return dataset_query_obj


def parse_stack_based_proc_section(root):
    """Parse the stack_based_processing section, for all the chains:
       AGB, FH, FD, TOMO FH"""

    stack_based_proc_item = root.find("stack_based_processing")

    if stack_based_proc_item:
        boundaries_whole = parse_boundaries_latlon_section(stack_based_proc_item)

        L1cProductList_Item = stack_based_proc_item.find("L1C_product_list")
        L1cProduct_Items = L1cProductList_Item.findall("L1C_product")

        stack_composition = {}

        ECEF_grid_file_names = {}
        kz_file_names = {}
        slant_range_distances_file_names = {}
        off_nadir_angle_file_names = {}
        slope_file_names = {}
        reference_height_file_names = {}
        calibration_screens_file_names = {}
        boundaries_per_stack = {}

        for L1cProduct_Item in L1cProduct_Items:

            stack_id = L1cProduct_Item.attrib["unique_stack_id"]

            # acquisition
            Acquisition_Items = L1cProduct_Item.findall("acquisition")
            Acquisitions_list = []
            for Acquisition_Item in Acquisition_Items:
                Acquisitions_list.append(Acquisition_Item.text)
            stack_composition[stack_id] = Acquisitions_list

            # geographic_boundaries (optional)
            boundaries_per_stack[stack_id] = parse_boundaries_latlon_section(L1cProduct_Item)

            # auxiliaries
            auxiliaries_item = L1cProduct_Item.find("auxiliaries")

            Geometry_Item = auxiliaries_item.find("geometry")

            ECEF_grid_file_names[stack_id] = Geometry_Item.find("ECEF_grid").text
            kz_file_names[stack_id] = Geometry_Item.find("vertical_wavenumber").text
            slant_range_distances_file_names[stack_id] = Geometry_Item.find("radar_distances").text
            off_nadir_angle_file_names[stack_id] = Geometry_Item.find("off_nadir_angles").text
            slope_file_names[stack_id] = Geometry_Item.find("slope").text
            reference_height_file_names[stack_id] = Geometry_Item.find("reference_height").text

            if auxiliaries_item.find("calibration_screens") != None:
                calibration_screens_file_names[stack_id] = auxiliaries_item.find("calibration_screens").text

        # "DEM" mandatory for all the chains:
        DEM_Item = stack_based_proc_item.find("DEM")
        dem_folder = DEM_Item.text

        # optionals, depending on the chain:
        average_covariance_folder = ""
        reference_agb_folder = ""
        system_decorrelation_fun_folder = ""
        forest_mask_catalogue_folder = ""
        forest_height_folder = ""

        # "FD"
        AVG_Item = stack_based_proc_item.find("average_covariance")
        if AVG_Item and AVG_Item.text:
            average_covariance_folder = AVG_Item.text

        # "AGB"
        ref_agb_item = stack_based_proc_item.find("reference_agb")
        if not ref_agb_item == None:
            reference_agb_folder = ref_agb_item.text
        fh_item = stack_based_proc_item.find("forest_height")
        if not fh_item == None:
            forest_height_folder = fh_item.text

        # FH
        sys_dec_fun_item = stack_based_proc_item.find("system_decorrelation_function")
        if not sys_dec_fun_item == None:
            system_decorrelation_fun_folder = sys_dec_fun_item.text
        else:
            system_decorrelation_fun_folder = ""

        # All but TOMO
        fm_item = stack_based_proc_item.find("forest_mask")
        if not fm_item == None:
            forest_mask_catalogue_folder = fm_item.text

        stack_based_proc_obj = stack_based_processing(
            stack_composition,
            ECEF_grid_file_names,
            kz_file_names,
            slant_range_distances_file_names,
            off_nadir_angle_file_names,
            slope_file_names,
            reference_height_file_names,
            calibration_screens_file_names,
            dem_folder,
            reference_agb_folder,
            forest_mask_catalogue_folder,
            system_decorrelation_fun_folder,
            average_covariance_folder,
            forest_height_folder,
            boundaries_whole,
            boundaries_per_stack,
        )

    else:
        stack_based_proc_obj = None

    return stack_based_proc_obj


def parse_core_proc_agb_section(root):

    core_proc_agb_item = root.find("core_processing_agb")

    if core_proc_agb_item:
        cal_item = core_proc_agb_item.find("lookup_table_cal")
        lut_cal = parse_lut_section(cal_item)

        fnf_item = core_proc_agb_item.find("lookup_table_fnf")
        lut_fnf = parse_lut_section(fnf_item)

        stacks_item = core_proc_agb_item.find("lookup_table_stacks")
        lut_stacks = parse_lut_section(stacks_item)

        core_proc_agb_obj = core_processing_agb(lut_cal, lut_fnf, lut_stacks)

    else:
        core_proc_agb_obj = None

    return core_proc_agb_obj


def parse_core_proc_fh_section(root, tomo_fh_flag=False):

    if tomo_fh_flag:
        core_proc_fh_item = root.find("core_processing_tomo_fh")
    else:
        core_proc_fh_item = root.find("core_processing_fh")

    if core_proc_fh_item:
        equi7_products_paths_item = core_proc_fh_item.find("equi7_products_paths")

        data_equi7_fnames = {}
        mask_equi7_fnames = {}
        for equi7_product_item in equi7_products_paths_item.findall("equi7_product"):
            stack_id = equi7_product_item.attrib["unique_stack_id"]
            data_equi7_fnames[stack_id] = []
            mask_equi7_fnames[stack_id] = []

            for equi7_tile_item in equi7_product_item.findall("equi7_tile"):

                data_equi7_fnames[stack_id].append(equi7_tile_item.find("path_data").text)
                mask_equi7_fnames[stack_id].append(equi7_tile_item.find("path_mask").text)

        core_proc_fh_obj = core_processing_fh(data_equi7_fnames, mask_equi7_fnames)

    else:
        core_proc_fh_obj = None

    return core_proc_fh_obj


def parse_core_proc_fd_section(root):
    """
       place older function to be developed yet 
    """

    core_proc_fd_item = root.find("core_processing_fd")

    if core_proc_fd_item:
        cycles_composition = {}  # dict which groupes together all the differente global cycles for s
        # tacks with same nominal geometry
        # cycles_composition[nominal_geometry_stack_id][global_cycle_number]
        cycles_composition_item = core_proc_fd_item.find("cycles_composition")

        for nominal_geometry_item in cycles_composition_item.findall("nominal_geometry"):

            nominal_geometry_stack_id = nominal_geometry_item.attrib["stack_id"]
            cycles_composition[nominal_geometry_stack_id] = {}

            for global_cycle_item in nominal_geometry_item.findall("global_cycle"):

                global_cycle_idx = int(global_cycle_item.attrib["cycle_number"])

                cycles_composition[nominal_geometry_stack_id][global_cycle_idx] = []

                for acquisition_item in global_cycle_item.findall("acquisition"):
                    cycles_composition[nominal_geometry_stack_id][global_cycle_idx].append(acquisition_item.text)

        core_proc_fd_obj = core_processing_fd(cycles_composition)
    else:

        core_proc_fd_obj = None

    return core_proc_fd_obj


def parse_core_proc_tomo_fh_section(root):

    # core processing section for FH and TOMO FH are the same
    core_proc_tomo_fh_obj = parse_core_proc_fh_section(root, tomo_fh_flag=True)

    return core_proc_tomo_fh_obj

def parse_core_proc_fnf_section(root):
    
    core_proc_fnf_item = root.find("core_processing_fnf")

    if core_proc_fnf_item:
        equi7_products_paths_item = core_proc_fnf_item.find("equi7_products_paths")
        
        data_equi7_fnames = {}
        mask_equi7_fnames = {}
        for equi7_product_item in equi7_products_paths_item.findall("equi7_product"):
            stack_id = equi7_product_item.attrib["unique_stack_id"]
            data_equi7_fnames[stack_id] = []
            mask_equi7_fnames[stack_id] = []
            
            for equi7_tile_item in equi7_product_item.findall("equi7_tile"):
                
                data_equi7_fnames[stack_id].append( equi7_tile_item.find("path_data").text )
                mask_equi7_fnames[stack_id].append( equi7_tile_item.find("path_mask").text )
            
        core_proc_fnf_obj = core_processing_fnf(data_equi7_fnames, mask_equi7_fnames)

    else:
        core_proc_fnf_obj = None

    return core_proc_fnf_obj


def write_configuration_file(conf_params_obj, conf_file_xml):
    """
    Write the configuration file:
    """

    root_item = Element("BioPALConfiguration")
    root_item.set("version", _XML_VERSION)

    write_gdal_section(root_item, conf_params_obj.gdal)
    write_proc_flags_section(root_item, conf_params_obj.processing_flags)
    write_ground_canc_section(root_item, conf_params_obj.ground_cancellation)
    write_est_agb_section(root_item, conf_params_obj.estimate_agb)
    write_est_fh_section(root_item, conf_params_obj.estimate_fh)
    write_est_tomo_fh_section(root_item, conf_params_obj.estimate_tomo_fh)
    write_change_det_fd_section(root_item, conf_params_obj.change_detection_fd)

    # write to file
    output_dir = os.path.dirname(conf_file_xml)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    ElementTree_indent(root_item)
    tree = ElementTree(root_item)
    tree.write(open(conf_file_xml, "w"), encoding="unicode")


def write_gdal_section(father_item, gdal_obj):

    if (gdal_obj.gdal_path is None or not gdal_obj.gdal_path) and (
        gdal_obj.gdal_environment_path is None or not gdal_obj.gdal_environment_path
    ):
        return

    gdal_item = SubElement(father_item, "gdal")
    gdal_path_item = SubElement(gdal_item, "gdal_path")
    gdal_environment_path_item = SubElement(gdal_item, "gdal_environment_path")

    if gdal_obj.gdal_path is None or not gdal_obj.gdal_path:
        gdal_path_item.text = ""
    else:
        gdal_path_item.text = gdal_obj.gdal_path

    if gdal_obj.gdal_environment_path is None or not gdal_obj.gdal_environment_path:
        gdal_environment_path_item.text = ""
    else:
        gdal_environment_path_item.text = gdal_obj.gdal_environment_path


def write_proc_flags_section(father_item, proc_flags_obj):

    proc_flags_item = SubElement(father_item, "processing_flags")

    enable_resampling = SubElement(proc_flags_item, "enable_resampling")
    enable_resampling.text = str(proc_flags_obj.enable_resampling)
    compute_geometry = SubElement(proc_flags_item, "compute_geometry")
    compute_geometry.text = str(proc_flags_obj.compute_geometry)
    apply_calibration_screen = SubElement(proc_flags_item, "apply_calibration_screen")
    apply_calibration_screen.text = str(proc_flags_obj.apply_calibration_screen)
    DEM_flattening = SubElement(proc_flags_item, "DEM_flattening")
    DEM_flattening.text = str(proc_flags_obj.DEM_flattening)
    multilook_heading_correction = SubElement(proc_flags_item, "multilook_heading_correction")
    multilook_heading_correction.text = str(proc_flags_obj.multilook_heading_correction)
    input_data_type = SubElement(proc_flags_item, "input_data_type")
    input_data_type.text = str(proc_flags_obj.input_data_type)
    save_breakpoints = SubElement(proc_flags_item, "save_breakpoints")
    save_breakpoints.text = str(proc_flags_obj.save_breakpoints)
    delete_temporary_files = SubElement(proc_flags_item, "delete_temporaryFiles")
    delete_temporary_files.text = str(proc_flags_obj.delete_temporary_files)


def write_ground_canc_section(father_item, ground_canc_obj):

    ground_cannc_item = SubElement(father_item, "ground_cancellation")
    MultiMaster = SubElement(ground_cannc_item, "multi_master")
    MultiMaster.text = str(ground_canc_obj.multi_master_flag)
    EnhancedForestHeight = SubElement(ground_cannc_item, "enhanced_forest_height")
    EnhancedForestHeight.text = str(ground_canc_obj.enhanced_forest_height)
    EnhancedForestHeight.set("unit", "m")
    ModelBasedEqualization = SubElement(ground_cannc_item, "model_based_equalization")
    ModelBasedEqualization.text = str(ground_canc_obj.equalization_flag)
    ModelBasedEqualization.set("always_off", "1")
    ModelBasedEqualization.set("always_on", "2")
    ModelBasedEqualization.set("on_if_two_acquisitions", "3")

    if (
        not ModelBasedEqualization.text == "1"
        and not ModelBasedEqualization.text == "2"
        and not ModelBasedEqualization.text == "3"
    ):
        error_str = 'Configuration flag "ModelBasedEqualization" value "{}" not valid, choose among "1", "2" or "3", where "1"->always OFF; "2"->always ON; "3"->ON only if two acquisitions are present'.format(
            ModelBasedEqualization.text
        )
        logging.error(error_str)
        raise ValueError(error_str)


def write_est_agb_section(father_item, est_agb_obj):

    est_agb_item = SubElement(father_item, "estimate_agb")

    product_resolution = SubElement(est_agb_item, "product_resolution")
    product_resolution.text = str(est_agb_obj.product_resolution)
    product_resolution.set("unit", "m")

    forest_class_observable_name = SubElement(est_agb_item, "forest_class_observable_name")
    forest_class_observable_name.text = str(est_agb_obj.forest_class_observable_name)

    transfer_function_name = SubElement(est_agb_item, "transfer_function_name")
    transfer_function_name.text = str(est_agb_obj.transfer_function_name)

    number_of_tests = SubElement(est_agb_item, "number_of_tests")
    number_of_tests.text = str(est_agb_obj.number_of_tests)

    fraction_of_roi_per_test = SubElement(est_agb_item, "fraction_of_roi_per_test")
    fraction_of_roi_per_test.text = str(est_agb_obj.fraction_of_roi_per_test)

    fraction_of_cal_per_test = SubElement(est_agb_item, "fraction_of_cal_per_test")
    fraction_of_cal_per_test.text = str(est_agb_obj.fraction_of_cal_per_test)

    add_variability_on_cal_data = SubElement(est_agb_item, "add_variability_on_cal_data")
    add_variability_on_cal_data.text = str(est_agb_obj.add_variability_on_cal_data)

    intermediate_ground_averaging = SubElement(est_agb_item, "intermediate_ground_averaging")
    intermediate_ground_averaging.text = str(est_agb_obj.intermediate_ground_averaging)

    distance_sampling_area = SubElement(est_agb_item, "distance_sampling_area")
    distance_sampling_area.text = str(est_agb_obj.distance_sampling_area)
    distance_sampling_area.set("unit", "m")

    parameter_block_size = SubElement(est_agb_item, "parameter_block_size")
    parameter_block_size.text = str(est_agb_obj.parameter_block_size)
    parameter_block_size.set("unit", "m")

    distance_parameter_block = SubElement(est_agb_item, "distance_parameter_block")
    distance_parameter_block.text = str(est_agb_obj.distance_parameter_block)
    distance_parameter_block.set("unit", "m")

    min_number_of_rois = SubElement(est_agb_item, "min_number_of_rois")
    min_number_of_rois.text = str(est_agb_obj.min_number_of_rois)

    min_number_of_rois_per_stack = SubElement(est_agb_item, "min_number_of_rois_per_stack")
    min_number_of_rois_per_stack.text = str(est_agb_obj.min_number_of_rois_per_stack)

    min_number_of_cals_per_test = SubElement(est_agb_item, "min_number_of_cals_per_test")
    min_number_of_cals_per_test.text = str(est_agb_obj.min_number_of_cals_per_test)

    min_number_of_rois_per_test = SubElement(est_agb_item, "min_number_of_rois_per_test")
    min_number_of_rois_per_test.text = str(est_agb_obj.min_number_of_rois_per_test)

    estimation_valid_values_limits = SubElement(est_agb_item, "estimation_valid_values_limits")
    min_item = SubElement(estimation_valid_values_limits, "min")
    min_item.text = str(est_agb_obj.estimation_valid_values_limits.min)
    min_item.set("unit", "t/ha")
    max_item = SubElement(estimation_valid_values_limits, "max")
    max_item.text = str(est_agb_obj.estimation_valid_values_limits.max)
    max_item.set("unit", "t/ha")

    write_residual_function_section(est_agb_item, est_agb_obj.residual_function)


def write_residual_function_section(Estimate_elem, res_fun_obj):

    residual_function = SubElement(Estimate_elem, "residual_function")

    formula_terms = SubElement(residual_function, "formula_terms")
    number_of_terms = len(res_fun_obj.formula_terms.string)
    term_struct = res_fun_obj.formula_terms
    for index in np.arange(number_of_terms):

        term = SubElement(formula_terms, "term")

        formula_name = SubElement(term, "name")
        formula_name.text = term_struct.name[index]

        formula_string = SubElement(term, "string")
        formula_string.text = term_struct.string[index]

        formula_weights = SubElement(term, "weight")

        formula_weight_step1 = SubElement(formula_weights, "step1")
        formula_weight_step1.text = str(term_struct.formula_weights.step1[index])
        formula_weight_step2 = SubElement(formula_weights, "step2")
        formula_weight_step2.text = str(term_struct.formula_weights.step2[index])
        formula_weight_step3 = SubElement(formula_weights, "step3")
        formula_weight_step3.text = str(term_struct.formula_weights.step3[index])

    formula_parameters = SubElement(residual_function, "formula_parameters")
    number_of_parameters = len(res_fun_obj.formula_parameters.name)
    par_struct = res_fun_obj.formula_parameters
    for index in np.arange(number_of_parameters):

        par = SubElement(formula_parameters, "par")

        name = SubElement(par, "name")
        name.text = par_struct.name[index]

        save_as_map = SubElement(par, "save_as_map")
        save_as_map.text = str(par_struct.save_as_map[index])

        transform = SubElement(par, "transform")
        transform.text = par_struct.transform[index]

        limits = SubElement(par, "limits")
        limits.set("unit", par_struct.limit_units[index])
        min_item = SubElement(limits, "min")
        min_item.text = str(par_struct.limits[index][0])
        max_item = SubElement(limits, "max")
        max_item.text = str(par_struct.limits[index][1])

        associated_observable = SubElement(par, "associated_observable_name")
        associated_observable.text = par_struct.associated_observable_name[index]

        parameter_variabilities = SubElement(par, "variability")
        samples_item = SubElement(parameter_variabilities, "samples")
        samples_item.text = str(par_struct.parameter_variabilities[index][0])
        forest_classes_item = SubElement(parameter_variabilities, "forest_classes")
        forest_classes_item.text = str(par_struct.parameter_variabilities[index][1])
        stacks_item = SubElement(parameter_variabilities, "stacks")
        stacks_item.text = str(par_struct.parameter_variabilities[index][2])
        global_cycles_item = SubElement(parameter_variabilities, "global_cycles")
        global_cycles_item.text = str(par_struct.parameter_variabilities[index][3])
        headings_item = SubElement(parameter_variabilities, "headings")
        headings_item.text = str(par_struct.parameter_variabilities[index][4])
        swaths_item = SubElement(parameter_variabilities, "swaths")
        swaths_item.text = str(par_struct.parameter_variabilities[index][5])
        subswaths_item = SubElement(parameter_variabilities, "subswaths")
        subswaths_item.text = str(par_struct.parameter_variabilities[index][6])
        azimuth_images_item = SubElement(parameter_variabilities, "azimuth_images")
        azimuth_images_item.text = str(par_struct.parameter_variabilities[index][7])

    formula_observables = SubElement(residual_function, "formula_observables")
    obs_struct = res_fun_obj.formula_observables
    # obs_struct.source composition:
    # obs_struct.source[index_obs][index_stack][index_file][index_layer]
    # and each layer has two fields: path and layer id
    # path = obs_struct.source[index_obs][index_stack][index_file][index_layer][0]
    # band_id = obs_struct.source[index_obs][index_stack][index_file][index_layer][0]
    number_of_observables = len(obs_struct.source_paths)
    for index_obs in np.arange(number_of_observables):
        source_curr = obs_struct.source_paths[index_obs]
        obs = SubElement(formula_observables, "obs")
        name = SubElement(obs, "name")
        name.text = obs_struct.name[index_obs]

        is_required = SubElement(obs, "is_required")
        is_required.text = str(obs_struct.is_required[index_obs])

        sources = SubElement(obs, "sources")
        sources.set("unit", str(obs_struct.source_unit[index_obs]))
        sources.set("resolution_m", str(obs_struct.source_resolution[index_obs]))
        number_of_stacks = len(source_curr)
        for index_stack in np.arange(number_of_stacks):
            stack_curr = obs_struct.source_paths[index_obs][index_stack]

            stack = SubElement(sources, "stack")
            number_of_files = len(stack_curr)
            for index_file in np.arange(number_of_files):
                file_curr = obs_struct.source_paths[index_obs][index_stack][index_file]

                file = SubElement(stack, "file")

                index_layer_text = 0
                index_band_id = 1
                layer = SubElement(file, "path")
                layer.text = file_curr[index_layer_text]
                layer.set("band", str(file_curr[index_band_id]))

        limits = SubElement(obs, "limits")
        if not obs_struct.limit_units[index_obs] == "":
            limits.set("unit", obs_struct.limit_units[index_obs])
        min_item = SubElement(limits, "min")
        min_item.text = str(obs_struct.limits[index_obs][0])
        max_item = SubElement(limits, "max")
        max_item.text = str(obs_struct.limits[index_obs][1])

        transform = SubElement(obs, "transform")
        transform.text = str(obs_struct.transform[index_obs])

        averaging_method = SubElement(obs, "averaging_method")
        averaging_method.text = str(obs_struct.averaging_method[index_obs])


def write_est_fh_section(father_item, est_fh_obj):
    """
       place older function to be developed yet , if required
    """
    pass


def write_est_tomo_fh_section(father_item, est_tomo_fh_obj):
    """
       place older function to be developed yet , if required
    """
    pass


def write_change_det_fd_section(father_item, change_det_fd_obj):
    """
       place older function to be developed yet , if required
    """
    pass


def parse_configuration_file(configuration_file_xml):

    tree = ET.parse(configuration_file_xml)
    root = tree.getroot()

    gdal_obj = parse_gdal_section(root)

    proc_flags_obj = parse_proc_flags_section(root)

    ground_canc_obj = parse_ground_canc_section(root)

    est_agb_obj = parse_est_agb_section(root)

    est_fh_obj = parse_est_fh_section(root.find("estimate_fh"))

    est_tomo_fh_obj = parse_est_tomo_fh_section(root.find("estimate_tomo_fh"))

    change_det_obj = parse_change_det_section(root.find("change_detection_fd"))
    
    est_fnf_obj = parse_estimate_fnf_section(root.find("estimate_fnf"))

    conf_params_obj = conf_params(
        gdal_obj, proc_flags_obj, ground_canc_obj, est_agb_obj, est_fh_obj,
        est_tomo_fh_obj, change_det_obj, est_fnf_obj
    )
    return conf_params_obj


def parse_gdal_section(father_item):

    gdal_item = father_item.find("gdal")
    if gdal_item:
        gdal_path = gdal_item.find("gdal_path").text
        gdal_environment_path = gdal_item.find("gdal_environment_path").text

        conf_gdal_obj = conf_gdal(gdal_path, gdal_environment_path)

    else:
        conf_gdal_obj = conf_gdal(None, None)

    return conf_gdal_obj


def parse_proc_flags_section(father_item):

    proc_flags_item = father_item.find("processing_flags")

    enable_resampling = bool_from_string(proc_flags_item.find("enable_resampling").text)
    compute_geometry = bool_from_string(proc_flags_item.find("compute_geometry").text)
    apply_calibration_screen = bool_from_string(proc_flags_item.find("apply_calibration_screen").text)
    DEM_flattening = bool_from_string(proc_flags_item.find("DEM_flattening").text)
    multilook_heading_correction = bool_from_string(proc_flags_item.find("multilook_heading_correction").text)
    input_data_type = proc_flags_item.find("input_data_type").text
    save_breakpoints = bool_from_string(proc_flags_item.find("save_breakpoints").text)
    delete_temporary_files = bool_from_string(proc_flags_item.find("delete_temporaryFiles").text)

    conf_flags_obj = conf_flags(
        enable_resampling,
        compute_geometry,
        apply_calibration_screen,
        DEM_flattening,
        multilook_heading_correction,
        input_data_type,
        save_breakpoints,
        delete_temporary_files,
    )

    return conf_flags_obj


def parse_ground_canc_section(father_item):

    ground_canc_item = father_item.find("ground_cancellation")
    if ground_canc_item:
        multi_master_flag = bool_from_string(ground_canc_item.find("multi_master").text)
        enhanced_forest_height = float(ground_canc_item.find("enhanced_forest_height").text)
        equalization_flag = ground_canc_item.find("model_based_equalization").text
        if not equalization_flag == "1" and not equalization_flag == "2" and not equalization_flag == "3":
            error_str = 'Configuration flag "ModelBasedEqualization" value "{}" not valid, choose among "1", "2" or "3", where "1"->always OFF; "2"->always ON; "3"->ON only if two acquisitions are present'.format(
                equalization_flag
            )
            logging.error(error_str)
            raise ValueError(error_str)

        ground_canc_obj = conf_ground_canc(multi_master_flag, enhanced_forest_height, equalization_flag)

    else:
        ground_canc_obj = None

    return ground_canc_obj


def parse_est_agb_section(father_item):

    estimate_agb_item = father_item.find("estimate_agb")

    if estimate_agb_item:

        product_resolution = float(estimate_agb_item.find("product_resolution").text)
        forest_class_observable_name = str(estimate_agb_item.find("forest_class_observable_name").text)
        transfer_function_name = str(estimate_agb_item.find("transfer_function_name").text)
        number_of_tests = int(estimate_agb_item.find("number_of_tests").text)
        fraction_of_roi_per_test = float(estimate_agb_item.find("fraction_of_roi_per_test").text)
        fraction_of_cal_per_test = float(estimate_agb_item.find("fraction_of_cal_per_test").text)
        add_variability_on_cal_data = bool_from_string(estimate_agb_item.find("add_variability_on_cal_data").text)
        intermediate_ground_averaging = float(estimate_agb_item.find("intermediate_ground_averaging").text)
        distance_sampling_area = float(estimate_agb_item.find("distance_sampling_area").text)
        parameter_block_size = float(estimate_agb_item.find("parameter_block_size").text)
        distance_parameter_block = float(estimate_agb_item.find("distance_parameter_block").text)
        min_number_of_rois = int(estimate_agb_item.find("min_number_of_rois").text)
        min_number_of_rois_per_stack = int(estimate_agb_item.find("min_number_of_rois_per_stack").text)
        min_number_of_cals_per_test = int(estimate_agb_item.find("min_number_of_cals_per_test").text)
        min_number_of_rois_per_test = int(estimate_agb_item.find("min_number_of_rois_per_test").text)
        estimation_min_value = float(estimate_agb_item.find("estimation_valid_values_limits").find("min").text)
        estimation_max_value = float(estimate_agb_item.find("estimation_valid_values_limits").find("max").text)
        estimation_valid_values_limits = min_max(estimation_min_value, estimation_max_value,)

        residual_function_obj = parse_agb_residual_function_section(estimate_agb_item)

        conf_agb_est_obj = conf_agb_est(
            product_resolution,
            forest_class_observable_name,
            transfer_function_name,
            number_of_tests,
            fraction_of_roi_per_test,
            fraction_of_cal_per_test,
            add_variability_on_cal_data,
            intermediate_ground_averaging,
            distance_sampling_area,
            parameter_block_size,
            distance_parameter_block,
            min_number_of_rois,
            min_number_of_rois_per_stack,
            min_number_of_cals_per_test,
            min_number_of_rois_per_test,
            estimation_valid_values_limits,
            residual_function_obj,
        )

    else:
        conf_agb_est_obj = None

    return conf_agb_est_obj


def parse_agb_residual_function_section(father_item):

    residual_function_item = father_item.find("residual_function")

    formula_Item = residual_function_item.find("formula_terms")
    formula_parameters_Item = residual_function_item.find("formula_parameters")
    formula_observables_Item = residual_function_item.find("formula_observables")

    formula_string = []
    formula_name = []
    formula_weight_step1 = []
    formula_weight_step2 = []
    formula_weight_step3 = []
    for term_item in formula_Item.findall("term"):
        formula_name.append(term_item.find("name").text)
        formula_string.append(term_item.find("string").text)
        formula_weight_step1.append(float(term_item.find("weight").find("step1").text))
        formula_weight_step2.append(float(term_item.find("weight").find("step2").text))
        formula_weight_step3.append(float(term_item.find("weight").find("step3").text))
    formula_weight_struct = formula_weights(formula_weight_step1, formula_weight_step2, formula_weight_step3)
    formula_terms_struct = formula_terms(formula_name, formula_string, formula_weight_struct,)

    name = []
    save_as_map = []
    transform = []
    limits = []
    units_par = []
    associated_observable = []
    variabilities = []
    for par_item in formula_parameters_Item.findall("par"):

        name.append(par_item.find("name").text)
        save_as_map.append(bool_from_string(par_item.find("save_as_map").text))
        transform.append(par_item.find("transform").text)
        min_limit = float(par_item.find("limits").find("min").text)
        max_limit = float(par_item.find("limits").find("max").text)
        units_par.append(par_item.find("limits").attrib["unit"])
        limits.append([min_limit, max_limit])
        associated_observable.append(par_item.find("associated_observable_name").text)
        variabilities.append(
            [
                bool_from_string(par_item.find("variability").find("samples").text),
                bool_from_string(par_item.find("variability").find("forest_classes").text),
                bool_from_string(par_item.find("variability").find("stacks").text),
                bool_from_string(par_item.find("variability").find("global_cycles").text),
                bool_from_string(par_item.find("variability").find("headings").text),
                bool_from_string(par_item.find("variability").find("swaths").text),
                bool_from_string(par_item.find("variability").find("subswaths").text),
                bool_from_string(par_item.find("variability").find("azimuth_images").text),
            ]
        )

    formula_parameters_struct = formula_parameters(
        name, save_as_map, transform, limits, units_par, associated_observable, variabilities,
    )

    name = []
    is_required = []
    source_paths = []  # contains one sub list for each observable:
    source_units = []
    source_resolutions = []
    # observables_sourcesobs_struct.source[index_obs][index_stack][index_file][index_layer]

    limits_units = []
    limits = []
    transform = []
    averaging_method = []
    for obs_item in formula_observables_Item.findall("obs"):

        name.append(obs_item.find("name").text)
        is_required.append(bool_from_string(obs_item.find("is_required").text))
        min_range = float(obs_item.find("limits").find("min").text)
        max_range = float(obs_item.find("limits").find("max").text)
        if obs_item.find("limits").attrib == {}:
            limits_units.append("")
        else:
            limits_units.append(obs_item.find("limits").attrib["unit"])
        limits.append([min_range, max_range])
        transform.append(obs_item.find("transform").text)
        averaging_method.append(obs_item.find("averaging_method").text)

        source_item = obs_item.find("sources")
        if obs_item.find("sources").attrib == {}:
            source_units.append("")
            source_resolutions.append(0)
        else:
            source_resolutions.append(float(obs_item.find("sources").attrib["resolution_m"]))
            source_units.append(obs_item.find("sources").attrib["unit"])

        curr_obs_stacks = []
        for stack_item in source_item.findall("stack"):

            curr_stack_files = []
            for file_item in stack_item.findall("file"):

                curr_file_layers = []
                for layer_item in file_item.findall("path"):

                    layer_path = layer_item.text
                    band_id = int(layer_item.attrib["band"])
                    curr_file_layers.append(layer_path)
                    curr_file_layers.append(band_id)

                curr_stack_files.append(curr_file_layers)
            curr_obs_stacks.append(curr_stack_files)
        source_paths.append(curr_obs_stacks)

    formula_observables_struct = formula_observables(
        name,
        is_required,
        source_paths,
        source_units,
        source_resolutions,
        limits_units,
        limits,
        transform,
        averaging_method,
    )

    residual_function_obj = conf_residual_function(
        formula_terms_struct, formula_parameters_struct, formula_observables_struct,
    )

    return residual_function_obj


def parse_est_fh_section(est_fh_item):

    if est_fh_item:

        product_resolution = float(est_fh_item.find("product_resolution").text)

        estimation_valid_values_limits = [
            float(est_fh_item.find("estimation_valid_values_limits").find("min").text),
            float(est_fh_item.find("estimation_valid_values_limits").find("max").text),
        ]

        spectral_shift_filtering = bool_from_string(est_fh_item.find("spectral_shift_filtering").text)

        kz_thresholds = [
            float(est_fh_item.find("KZ_thresholds").find("min").text),
            float(est_fh_item.find("KZ_thresholds").find("max").text),
        ]

        median_factor_temp = float(est_fh_item.find("median_factor").text)
        median_factor = int(median_factor_temp)
        if median_factor_temp - median_factor > 0:
            logging.warning(
                "FH configuration file MedianFactor should be integer, the value that will be used is: {}".format(
                    median_factor
                )
            )

        model_parameters_Item = est_fh_item.find("model_parameters")
        maximum_height = int(model_parameters_Item.find("maximum_height").text)
        number_of_extinction_value = int(model_parameters_Item.find("number_of_extinction_value").text)
        number_of_ground_volume_ratio_value = int(
            model_parameters_Item.find("number_of_ground_volume_ratio_value").text
        )
        number_of_temporal_decorrelation_value = int(
            model_parameters_Item.find("number_of_temporal_decorrelation_value").text
        )
        if (
            maximum_height < 2
            or number_of_extinction_value < 2
            or number_of_ground_volume_ratio_value < 2
            or number_of_temporal_decorrelation_value < 2
        ):
            error_message = 'Following FH configurations should be all values GREATER than "1": "MaximumHeight", "NumberOfExtinctionValue", "NumberOfGroundVolumeRatioValue" and "NumberOfTemporalDecorrelationValue"'
            logging.error(error_message)
            raise ValueError(error_message)
        model_parameters = conf_fh_model_params(
            estimation_valid_values_limits,
            maximum_height,
            number_of_extinction_value,
            number_of_ground_volume_ratio_value,
            number_of_temporal_decorrelation_value,
        )

        conf_fh_est_obj = conf_fh_est(
            spectral_shift_filtering, product_resolution, kz_thresholds, model_parameters, median_factor,
        )

    else:
        conf_fh_est_obj = None

    return conf_fh_est_obj


def parse_change_det_section(change_detect_item):

    if change_detect_item:

        product_resolution = float(change_detect_item.find("product_resolution").text)
        confidence_level = float(change_detect_item.find("confidence_level").text)

        conf_fd_est_obj = conf_fd_est(product_resolution, confidence_level)

    else:
        conf_fd_est_obj = None

    return conf_fd_est_obj


def parse_estimate_fnf_section(estimate_fnf_item):

    if estimate_fnf_item:

        product_resolution = float(estimate_fnf_item.find("product_resolution").text)
        logreg_coeffs_fname = (estimate_fnf_item.find("logreg_coeffs_fname").text)

        conf_fnf_est_obj = conf_fnf_est(product_resolution, logreg_coeffs_fname)

    else:
        conf_fnf_est_obj = None

    return conf_fnf_est_obj


def parse_est_tomo_fh_section(tomo_fh_item):

    if tomo_fh_item:

        product_resolution = float(tomo_fh_item.find("product_resolution").text)

        VerticalRange_Item = tomo_fh_item.find("vertical_range")
        maximum_height = float(VerticalRange_Item.find("maximum_height").text)
        minimum_height = float(VerticalRange_Item.find("minimum_height").text)
        sampling = float(VerticalRange_Item.find("sampling").text)
        vertical_range = vertical_range_params(maximum_height, minimum_height, sampling)

        estimation_valid_values_limits = [
            float(tomo_fh_item.find("estimation_valid_values_limits").find("min").text),
            float(tomo_fh_item.find("estimation_valid_values_limits").find("max").text),
        ]

        enable_super_resolution = bool_from_string(tomo_fh_item.find("enable_super_resolution").text)
        regularization_noise_factor = float(tomo_fh_item.find("regularization_noise_factor").text)
        power_threshold = float(tomo_fh_item.find("power_threshold").text)
        median_factor = int(tomo_fh_item.find("median_factor").text)

        conf_tomo_fh_est_obj = conf_tomo_fh_est(
            product_resolution,
            vertical_range,
            estimation_valid_values_limits,
            enable_super_resolution,
            regularization_noise_factor,
            power_threshold,
            median_factor,
        )

    else:
        conf_tomo_fh_est_obj = None

    return conf_tomo_fh_est_obj


############## core functions
def ElementTree_indent(elem, level=0):
    i = "\n" + level * "  "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for elem in elem:
            ElementTree_indent(elem, level + 1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i


def bool_from_string(in_string):

    in_string = in_string[0].upper() + in_string[1:].lower()

    if "True" in in_string or "False" in in_string:
        bool_value = ast.literal_eval(in_string)

    elif "0" in in_string or "1" in in_string:
        bool_value = bool(ast.literal_eval(in_string))
    else:
        error_message = "Main input file: the L2product AGB, FH, FD and TOMO values should be boolean: True, False, 0 or 1 are supported"
        logging.error(error_message)
        raise ValueError(error_message)

    return bool_value


################ Boundaries sections core parsers and readers ################
def write_boundaries_latlon_section(father_item, geographic_boundaries_struct):

    geographic_boundaries_Item = SubElement(father_item, "geographic_boundaries")
    geographic_boundaries_Item.set("unit", "deg")

    latitude_Item = SubElement(geographic_boundaries_Item, "latitude")
    min_latitude_Item = SubElement(latitude_Item, "min")
    min_latitude_Item.text = str(geographic_boundaries_struct.lat_min)
    max_latitude_Item = SubElement(latitude_Item, "max")
    max_latitude_Item.text = str(geographic_boundaries_struct.lat_max)

    longitude_Item = SubElement(geographic_boundaries_Item, "longitude")
    min_longitude_Item = SubElement(longitude_Item, "min")
    min_longitude_Item.text = str(geographic_boundaries_struct.lon_min)
    max_longitude_Item = SubElement(longitude_Item, "max")
    max_longitude_Item.text = str(geographic_boundaries_struct.lon_max)

    return father_item


def write_boundaries_eastnorth_section(father_item, boundaries):

    geographic_boundaries_Item = SubElement(father_item, "geographic_boundaries")
    geographic_boundaries_Item.set("unit", "m")

    east_Item = SubElement(geographic_boundaries_Item, "east")
    min_east_Item = SubElement(east_Item, "min")
    min_east_Item.text = str(boundaries[0])
    max_east_Item = SubElement(east_Item, "max")
    max_east_Item.text = str(boundaries[1])

    north_Item = SubElement(geographic_boundaries_Item, "north")
    min_north_Item = SubElement(north_Item, "min")
    min_north_Item.text = str(boundaries[2])
    max_north_Item = SubElement(north_Item, "max")
    max_north_Item.text = str(boundaries[3])

    if len(boundaries) == 5:
        flag_cal = SubElement(father_item, "flag_cal_format")
        flag_cal.set("formats", "1= RASTER ; 0=GeoJSON")
        flag_cal.text = str(int(boundaries[4]))


def parse_boundaries_latlon_section(fatherItem):
    """
        This function is shared by all the parser containing boundaries: 
            this is an optional field, depending on the processing chain, so a 
            check is done here, and if absent, None is returned
    """
    boundaries_out = geographic_boundaries_latlon(None, None, None, None)

    boundaries_Item = fatherItem.find("geographic_boundaries")
    if boundaries_Item:
        boundaries_out.lon_min = float(boundaries_Item.find("longitude").find("min").text)
        boundaries_out.lon_max = float(boundaries_Item.find("longitude").find("max").text)
        boundaries_out.lat_min = float(boundaries_Item.find("latitude").find("min").text)
        boundaries_out.lat_max = float(boundaries_Item.find("latitude").find("max").text)

    return boundaries_out


def parse_boundaries_eastnorth_section(fatherItem):
    """
        This function is shared by all the parser containing boundaries
    """
    boundaries_out = geographic_boundaries_northeast(None, None, None, None)

    boundaries_Item = fatherItem.find("geographic_boundaries")
    if boundaries_Item:
        boundaries_out.east_min = float(boundaries_Item.find("east").find("min").text)
        boundaries_out.east_max = float(boundaries_Item.find("east").find("max").text)
        boundaries_out.north_min = float(boundaries_Item.find("north").find("min").text)
        boundaries_out.north_max = float(boundaries_Item.find("north").find("max").text)

    return boundaries_out


def parse_polygon_boundaries_section(fatherItem):

    point_list = []
    GeographicBoundary_Item = fatherItem.find("geographic_boundaries_polygon")
    Points = GeographicBoundary_Item.findall("point")
    if len(Points) < 3:
        error_message = "Main input file: at least #3 Geograplic Boundary points should be specified"
        logging.error(error_message)
        raise ValueError(error_message)

    for Point_Item in Points:

        Latitude_str = Point_Item.find("latitude").text
        Longitude_str = Point_Item.find("longitude").text
        Lat = float(Latitude_str)
        Lon = float(Longitude_str)
        if Lat < -90 or Lat > 90:
            error_message = "Main input file: Latitude should be a number from -90 to +90 [deg] "
            logging.error(error_message)
            raise ValueError(error_message)
        if Lon < -180 or Lon > 180:
            error_message = "Main input file: Longitude should be a number from -180 to +180 [deg] "
            logging.error(error_message)
            raise ValueError(error_message)

        point_dict = {"Latitude": Lat, "Longitude": Lon}

        point_list.append(point_dict)
        del point_dict

    geographic_boundary_struct = geographic_boundary(point_list)

    return geographic_boundary_struct


############# core processing core lut parser and writer
def write_lut_section(fatherItem, lut, lut_type):
    """Write lookup table for lut_cal, lut_fnf and lut_stacks in xml format"""

    if not (lut_type == "cal" or lut_type == "fnf" or lut_type == "stacks"):
        errormessage = "type {} not recognized".format(lut_type)
        raise ValueError(errormessage)

    lut_item = SubElement(fatherItem, "lookup_table_" + lut_type)

    num_of_items = len(lut.paths)
    for item_idx in np.arange(num_of_items):
        item_item = SubElement(lut_item, "item")
        write_boundaries_eastnorth_section(item_item, lut.boundaries[item_idx, :])

        path_item = SubElement(item_item, "path")
        path_item.text = str(lut.paths[item_idx])

        if not lut.progressive is None:
            progressive_item = SubElement(item_item, "progressive_stacks_indexes")

            stack_progressive_index = SubElement(progressive_item, "stack_progressive_index")
            stack_progressive_index.text = str(int(lut.progressive[item_idx][0]))
            global_cycle_index = SubElement(progressive_item, "global_cycle_index")
            global_cycle_index.text = str(int(lut.progressive[item_idx][1]))
            heading = SubElement(progressive_item, "heading")
            heading.set("unit", "deg")
            heading.text = str(lut.progressive[item_idx][2])
            range_swath_index = SubElement(progressive_item, "range_swath_index")
            range_swath_index.text = str(int(lut.progressive[item_idx][3]))
            range_sub_swath_index = SubElement(progressive_item, "range_sub_swath_index")
            range_sub_swath_index.text = str(int(lut.progressive[item_idx][4]))
            azimuth_swath_index = SubElement(progressive_item, "azimuth_swath_index")
            azimuth_swath_index.text = str(int(lut.progressive[item_idx][5]))


def parse_lut_section(fatherItem):

    first_item = fatherItem.find("item")

    lut_type = "fnf"
    if not first_item.find("progressive_stacks_indexes") == None:
        lut_type = "stack"
    if not first_item.find("flag_cal_format") == None:
        lut_type = "cal"

    # initializations
    number_of_stacks = len(fatherItem.findall("item"))
    lut_paths = []
    if lut_type == "cal":
        lut_boundaries = np.zeros((number_of_stacks, 5))  # additional element (flag cal)
        lut_progressive = None
    if lut_type == "fnf":
        lut_boundaries = np.zeros((number_of_stacks, 4))
        lut_progressive = None
    if lut_type == "stack":
        lut_boundaries = np.zeros((number_of_stacks, 4))
        lut_progressive = np.zeros((number_of_stacks, 6))
    for stack_idx, lut_item in enumerate(fatherItem.findall("item")):

        # fill lut_paths
        lut_paths.append(str(lut_item.find("path").text))

        # fill lut_boundaries
        boundaries_curr = parse_boundaries_eastnorth_section(lut_item)

        if lut_type == "cal":

            flag_cal_format = lut_item.find("flag_cal_format").text
            lut_boundaries[stack_idx, :] = [
                float(boundaries_curr.east_min),
                float(boundaries_curr.east_max),
                float(boundaries_curr.north_min),
                float(boundaries_curr.north_max),
                int(flag_cal_format),
            ]
        else:
            lut_boundaries[stack_idx, :] = [
                float(boundaries_curr.east_min),
                float(boundaries_curr.east_max),
                float(boundaries_curr.north_min),
                float(boundaries_curr.north_max),
            ]

        # fill lut_progressive
        if lut_type == "stack":
            progressive_stack_item = lut_item.find("progressive_stacks_indexes")
            lut_progressive[stack_idx, :] = [
                int(float(progressive_stack_item.find("stack_progressive_index").text)),
                int(float(progressive_stack_item.find("global_cycle_index").text)),
                float(progressive_stack_item.find("heading").text),
                int(float(progressive_stack_item.find("range_swath_index").text)),
                int(float(progressive_stack_item.find("range_sub_swath_index").text)),
                int(float(progressive_stack_item.find("azimuth_swath_index").text)),
            ]

    lut_parsed = lut(lut_paths, lut_boundaries, lut_progressive)

    return lut_parsed
