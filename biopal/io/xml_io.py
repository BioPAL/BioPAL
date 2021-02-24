import os
import ast
import copy
import logging
import numpy as np
from lxml import etree as xml_tools
import xml.etree.ElementTree as ET
from xml.etree.ElementTree import ElementTree, Element, SubElement
from collections import namedtuple
from namedlist import namedlist
from arepytools.timing.precisedatetime import PreciseDateTime


###############################################################################
# structures (namedtuples) for read and write inputs and configuration files ##
_XML_VERSION = "1.0"

geographic_boundaries = namedlist(
    "geographic_boundaries",
    "lon_min \
                     lon_max \
                     lat_min \
                     lat_max",
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
 lines_start_utc",
)

# Main Input parameters structure (namedtuple)
main_input_params = namedtuple(
    "main_input_params",
    "L1c_repository \
 L1c_aux_data_repository \
 proc_flags \
 L1c_date \
 geographic_boundary \
 output_folder \
 geographic_grid_sampling \
 stacks_to_find",
)

# main_input_params "proc_flags" sub-fields:
proc_flags = namedlist(
    "proc_flags_params",
    "AGB \
                         FH \
                         TOMO_FH \
                         FD \
                         TOMO",
)

# main_input_params "L1c_date" sub-fields:
L1c_date = namedtuple(
    "L1c_date_params",
    "start \
                       stop",
)

# main_input_params "L1c_date" sub-fields:
geographic_boundary = namedtuple("geographic_boundary_params", "point")

# Input parameters structure (namedtuple)
input_params = namedtuple(
    "input_params",
    "chain_id \
 L1c_repository \
 stack_composition \
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
 output_folder \
 geographic_grid_sampling \
 system_decorrelation_fun_folder \
 average_covariance_folder \
 forest_height_folder",
)

# Configuration parameters structure (namedtuple)
configuration_params = namedtuple(
    "configuration_params",
    "chain_id \
 interpolate_stack \
 ground_cancellation \
 vertical_range \
 AGB \
 FH \
 FD \
 TOMO_FH \
 TOMO \
 enable_resampling \
 compute_geometry \
 apply_calibration_screen \
 DEM_flattening \
 multilook_heading_correction \
 save_breakpoints \
 delete_temporary_files",
)

# define each sub-structure of the configuration parameters, with sub-namedtuples
# configuration_params "interpolate_stack"  sub-fields:
interpolate_params = namedtuple("interpolate_params", "regular_baseline_grid_flag")
# configuration_params "ground_cancellation"  sub-fields:
ground_canc_params = namedtuple("ground_canc_params", "multi_master_flag enhanced_forest_height equalization_flag")

# configuration_params "agb" sub-fields:
agb_est_params = namedtuple(
    "agb_est_params",
    "residual_function \
                                 number_of_tests \
                                 forest_class_observable_name \
                                 transfer_function_name \
                                 fraction_of_roi_per_test \
                                 fraction_of_cal_per_test \
                                 add_variability_on_cal_data \
                                 intermediate_ground_averaging \
                                 product_resolution \
                                 distance_sampling_area \
                                 parameter_block_size \
                                 distance_parameter_block \
                                 min_number_of_rois \
                                 min_number_of_rois_per_stack \
                                 min_number_of_cals_per_test \
                                 min_number_of_rois_per_test \
                                 estimation_valid_values_limits",
)
# agb "residual_function" sub-fields:
residual_function = namedtuple("residual_function", "formula_terms formula_parameters formula_observables")
# agb "formula_terms" sub-fields:
formula_terms = namedtuple(
    "formula_terms",
    "name \
     string \
     formula_weights")
     
formula_weights = namedtuple(
    "formula_weights",
    "step1 \
     step2 \
     step3")
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
FH_est_params = namedtuple(
    "FH_est_params",
    "spectral_shift_filtering \
                                product_resolution \
                                kz_thresholds \
                                model_parameters \
                                median_factor",
)
# FH "model_parameters" sub-fields:
FH_model_params = namedtuple(
    "FH_model_params",
    "estimation_valid_values_limits maximum_height number_of_extinction_value number_of_ground_volume_ratio_value number_of_temporal_decorrelation_value",
)

# configuration_params "FD" sub-fields:
FD_est_params = namedtuple(
    "FD_est_params",
    "product_resolution \
                                confidence_level",
)

# configuration_params "TOMO FH" sub-fields:
TOMO_FH_est_params = namedtuple(
    "TOMO_FH_est_params",
    "estimation_valid_values_limits \
                                 product_resolution \
                                 enable_super_resolution \
                                 regularization_noise_factor \
                                 power_threshold \
                                 median_factor",
)

# configuration_params "TOMO" sub-fields:
TOMO_est_params = namedtuple("TOMO_est_params", "product_resolution")

# TOMO and TOMO_FH "vertical range" sub-fields:
vertical_range_params = namedtuple("vertical_range_params", "maximum_height minimum_height sampling")
###############################################################################


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


def write_boundaries_files(geographic_boundaries_whole, geographic_boundaries_per_stack, output_folder):
    """Write the two xmls containing boundaries"""

    # 1 of 2) boundaries file
    boundaries_fname = os.path.join(output_folder, "geographic_boundaries.xml")

    boundaries_Item = Element("boundaries")
    boundaries_Item.set("version", _XML_VERSION)

    boundaries_Item = write_boundaries_latlon_core(boundaries_Item, geographic_boundaries_whole)

    ElementTree_indent(boundaries_Item)
    tree = ElementTree(boundaries_Item)
    tree.write(open(boundaries_fname, "w"), encoding="unicode")

    # 2 of 2) boundaries per stack file
    boundaries_per_stack_fname = os.path.join(output_folder, "geographic_boundaries_per_stack.xml")

    boundaries_per_stack_Item = Element("boundaries_per_stack")
    boundaries_per_stack_Item.set("version", _XML_VERSION)

    for stack_id, geographic_boundaries_curr in geographic_boundaries_per_stack.items():
        stack_Item = SubElement(boundaries_per_stack_Item, "stack")
        stack_Item.set("stack_id", stack_id)
        stack_Item = write_boundaries_latlon_core(stack_Item, geographic_boundaries_curr)

    ElementTree_indent(boundaries_per_stack_Item)
    tree = ElementTree(boundaries_per_stack_Item)
    tree.write(open(boundaries_per_stack_fname, "w"), encoding="unicode")

    return boundaries_fname, boundaries_per_stack_fname


def write_boundaries_latlon_core(boundaries_Item, geographic_boundaries_struct):
    latitude_Item = SubElement(boundaries_Item, "latitude")
    latitude_Item.set("units", "Deg")
    max_latitude_Item = SubElement(latitude_Item, "max")
    max_latitude_Item.text = str(geographic_boundaries_struct.lat_max)
    min_latitude_Item = SubElement(latitude_Item, "min")
    min_latitude_Item.text = str(geographic_boundaries_struct.lat_min)

    longitude_Item = SubElement(boundaries_Item, "longitude")
    longitude_Item.set("units", "Deg")
    max_longitude_Item = SubElement(longitude_Item, "max")
    max_longitude_Item.text = str(geographic_boundaries_struct.lon_max)
    min_longitude_Item = SubElement(longitude_Item, "min")
    min_longitude_Item.text = str(geographic_boundaries_struct.lon_min)

    return boundaries_Item


def parse_boundaries_files(boundaries_fname):
    """Parse the two xmls containing boundaries (input can be the file name of 
    geographic_boundaries or geographic_boundaries_per_stack file name)"""

    tree = ET.parse(boundaries_fname)
    root = tree.getroot()

    if root.tag == "boundaries":
        geographic_boundaries_out = geographic_boundaries(None, None, None, None)
        geographic_boundaries_out.lon_min = float(root.find("longitude").find("min").text)
        geographic_boundaries_out.lon_max = float(root.find("longitude").find("max").text)
        geographic_boundaries_out.lat_min = float(root.find("latitude").find("min").text)
        geographic_boundaries_out.lat_max = float(root.find("latitude").find("max").text)

    elif root.tag == "boundaries_per_stack":
        geographic_boundaries_out = {}
        for stack_Item in root.findall("stack"):
            stack_id = stack_Item.attrib["stack_id"]

            geographic_boundaries_out[stack_id] = geographic_boundaries(None, None, None, None)
            geographic_boundaries_out[stack_id].lon_min = float(stack_Item.find("longitude").find("min").text)
            geographic_boundaries_out[stack_id].lon_max = float(stack_Item.find("longitude").find("max").text)
            geographic_boundaries_out[stack_id].lat_min = float(stack_Item.find("latitude").find("min").text)
            geographic_boundaries_out[stack_id].lat_max = float(stack_Item.find("latitude").find("max").text)

    return geographic_boundaries_out


def write_lut_files(lut, lut_type, output_folder):
    """Write loolup table for lut_cal, lut_fnf and lut_stacks in xml format"""

    if not (lut_type == "cal" or lut_type == "fnf" or lut_type == "stacks"):
        errormessage = "type {} not recognized".format(lut_type)
        raise ValueError(errormessage)

    lut_fname = os.path.join(output_folder, "lut_" + lut_type + ".xml")

    lut_Item = Element("lookpu_table_" + lut_type)
    lut_Item.set("version", _XML_VERSION)

    boundaries_elem = SubElement(lut_Item, "boundaries")
    boundaries_elem.set("units", "m")
    for boundaries in lut.boundaries:
        coordinates_elem = SubElement(boundaries_elem, "boundary_coordinates")

        east_min = SubElement(coordinates_elem, "east_min")
        east_min.text = str(boundaries[0])
        east_max = SubElement(coordinates_elem, "east_max")
        east_max.text = str(boundaries[1])
        north_min = SubElement(coordinates_elem, "north_min")
        north_min.text = str(boundaries[2])
        north_max = SubElement(coordinates_elem, "north_max")
        north_max.text = str(boundaries[3])

        if len(boundaries) == 5:
            flag_cal = SubElement(coordinates_elem, "flag_cal_format")
            flag_cal.set("formats", "1= RASTER ; 0=GeoJSON")
            flag_cal.text = str(int(boundaries[4]))

    paths_elem = SubElement(lut_Item, "paths")
    for path_curr in lut.paths:
        path_elem = SubElement(paths_elem, "path")
        path_elem.text = str(path_curr)

    if not lut.progressive is None:
        progressive_Item = SubElement(lut_Item, "progressive_stacks_indexes")
        for progressive_stack in lut.progressive:
            progressive_stack_elem = SubElement(progressive_Item, "progressive_stack")
            stack_progressive_index = SubElement(progressive_stack_elem, "stack_progressive_index")
            stack_progressive_index.text = str(int(progressive_stack[0]))
            global_cycle_index = SubElement(progressive_stack_elem, "global_cycle_index")
            global_cycle_index.text = str(int(progressive_stack[1]))
            heading = SubElement(progressive_stack_elem, "heading")
            heading.set("units", "Deg")
            heading.text = str(progressive_stack[2])
            range_swath_index = SubElement(progressive_stack_elem, "range_swath_index")
            range_swath_index.text = str(int(progressive_stack[3]))
            range_sub_swath_index = SubElement(progressive_stack_elem, "range_sub_swath_index")
            range_sub_swath_index.text = str(int(progressive_stack[4]))
            azimuth_swath_index = SubElement(progressive_stack_elem, "azimuth_swath_index")
            azimuth_swath_index.text = str(int(progressive_stack[5]))

    ElementTree_indent(lut_Item)
    tree = ElementTree(lut_Item)

    tree.write(open(lut_fname, "w"), encoding="unicode")


def parse_lut_files(lut_fname):

    tree = ET.parse(lut_fname)
    root = tree.getroot()

    boundaries_Item = root.find("boundaries")
    number_of_stacks = len(boundaries_Item.findall("boundary_coordinates"))

    if boundaries_Item.find("boundary_coordinates").find("flag_cal_format") == None:
        lut_type = "not_cal"
    else:
        lut_type = "cal"

    if lut_type == "cal":
        lut_boundaries = np.zeros((number_of_stacks, 5))  # additional element (flag cal)
    else:
        lut_boundaries = np.zeros((number_of_stacks, 4))

    for stack_idx, coord_item in enumerate(boundaries_Item.findall("boundary_coordinates")):
        if lut_type == "cal":
            lut_boundaries[stack_idx, :] = [
                float(coord_item.find("east_min").text),
                float(coord_item.find("east_max").text),
                float(coord_item.find("north_min").text),
                float(coord_item.find("north_max").text),
                int(coord_item.find("flag_cal_format").text),
            ]
        else:
            lut_boundaries[stack_idx, :] = [
                float(coord_item.find("east_min").text),
                float(coord_item.find("east_max").text),
                float(coord_item.find("north_min").text),
                float(coord_item.find("north_max").text),
            ]

    paths_Item = root.find("paths")
    lut_paths = []
    for path_elem in paths_Item.findall("path"):
        lut_paths.append(str(path_elem.text))

    progressive_Item = root.find("progressive_stacks_indexes")
    if progressive_Item:
        lut_progressive = np.zeros((number_of_stacks, 6))
        for stack_idx, progressive_stack_item in enumerate(progressive_Item.findall("progressive_stack")):
            lut_progressive[stack_idx, :] = [
                int(float(progressive_stack_item.find("stack_progressive_index").text)),
                int(float(progressive_stack_item.find("global_cycle_index").text)),
                float(progressive_stack_item.find("heading").text),
                int(float(progressive_stack_item.find("range_swath_index").text)),
                int(float(progressive_stack_item.find("range_sub_swath_index").text)),
                int(float(progressive_stack_item.find("azimuth_swath_index").text)),
            ]
    else:
        lut_progressive = None

    return lut_paths, lut_boundaries, lut_progressive


def write_chains_input_file(input_params, input_file_xml):
    """Write the input file of the inner chains"""

    chain_id = input_params.chain_id
    InputsL2_Item = Element("InputsL2" + chain_id)
    InputsL2_Item.set("version", _XML_VERSION)

    L1cRepository = SubElement(InputsL2_Item, "L1cRepository")
    L1cRepository.text = input_params.L1c_repository

    # L1cProductList
    L1cProductList = SubElement(InputsL2_Item, "L1cProductList")

    for stack_id, stack_composition_list in input_params.stack_composition.items():

        L1cProduct = SubElement(L1cProductList, "L1cProduct")
        L1cProduct.set("unique_stack_id", stack_id)

        for pf_name in stack_composition_list:

            Acquisition = SubElement(L1cProduct, "Acquisition")
            Acquisition.text = pf_name

    # AuxiliaryProductList
    AuxiliaryProductList = SubElement(InputsL2_Item, "AuxiliaryProductList")

    if len(input_params.ECEF_grid_file_names):
        for stack_id, stack_composition_list in input_params.stack_composition.items():

            L1cStack = SubElement(AuxiliaryProductList, "L1cStack")
            L1cStack.set("unique_stack_id", stack_id)

            Geometry = SubElement(L1cStack, "Geometry")

            ECEFGrid = SubElement(Geometry, "ECEFGrid")
            ECEFGrid.text = input_params.ECEF_grid_file_names[stack_id]

            kz = SubElement(Geometry, "kz")
            kz.text = input_params.kz_file_names[stack_id]

            RadarDistances = SubElement(Geometry, "RadarDistances")
            RadarDistances.text = input_params.slant_range_distances_file_names[stack_id]

            OffNadirAngle = SubElement(Geometry, "OffNadirAngles")
            OffNadirAngle.text = input_params.off_nadir_angle_file_names[stack_id]

            Slope = SubElement(Geometry, "Slope")
            Slope.text = input_params.slope_file_names[stack_id]

            ReferenceHeight = SubElement(Geometry, "ReferenceHeight")
            ReferenceHeight.text = input_params.reference_height_file_names[stack_id]

            CalibrationScreens = SubElement(L1cStack, "CalibrationScreens")
            CalibrationScreens.text = input_params.calibration_screens_file_names[stack_id]

    if chain_id == "FD":
        AverageCovariance = SubElement(AuxiliaryProductList, "AverageCovarianceFolder")
        AverageCovariance.text = input_params.average_covariance_folder

    DEM = SubElement(AuxiliaryProductList, "DEM")
    DEM.text = input_params.dem_folder

    if chain_id == "AGB":
        ReferenceAGB = SubElement(AuxiliaryProductList, "ReferenceAGB")
        ReferenceAGB.text = input_params.reference_agb_folder

        ForestHeight = SubElement(AuxiliaryProductList, "ForestHeight")
        ForestHeight.text = input_params.forest_height_folder

    if chain_id == "FH":
        SystemDecorrelationFunction = SubElement(AuxiliaryProductList, "SystemDecorrelationFunction")
        SystemDecorrelationFunction.text = input_params.system_decorrelation_fun_folder

    if not "TOMO" in chain_id:
        ForestMask = SubElement(AuxiliaryProductList, "ForestMask")
        ForestMask.text = input_params.forest_mask_catalogue_folder

    # OutputFolder
    OutputFolder = SubElement(InputsL2_Item, "OutputFolder")
    OutputFolder.text = input_params.output_folder

    # GeographicGridSampling
    GeographicGridSampling = SubElement(InputsL2_Item, "GeographicGridSampling")
    GeographicGridSampling.text = str(input_params.geographic_grid_sampling)
    GeographicGridSampling.set("unit", "m")

    # write to file
    ElementTree_indent(InputsL2_Item)
    tree = ElementTree(InputsL2_Item)
    tree.write(open(input_file_xml, "w"), encoding="unicode")


def write_chains_configuration_file(configuration_params, configuration_file_xml):
    """Write the configuration file of the inner chains FD FH TOMO and TOMO FH, not AGB"""

    chain_id = configuration_params.chain_id
    if chain_id == "AGB":
        error_str = "To write the AGB configuratin files, use its own specific functions"
        logging.error(error_str)
        raise ValueError(error_str)

    ConfigurationL2_Item = Element("ConfigurationL2" + chain_id)
    ConfigurationL2_Item.set("version", _XML_VERSION)

    if chain_id == "TOMO":
        InterpolateStack = SubElement(ConfigurationL2_Item, "InterpolateStack")
        RegularBaselineGrid = SubElement(InterpolateStack, "RegularBaselineGrid")
        RegularBaselineGrid.text = str(configuration_params.interpolate_stack.regular_baseline_grid_flag)

    # Ground Cancellation
    if chain_id == "FD":

        ground_canc_params = configuration_params.ground_cancellation
        ConfigurationL2_Item = write_ground_cancellation_core(ConfigurationL2_Item, ground_canc_params)

    params_curr = getattr(configuration_params, chain_id)

    if chain_id == "FH":

        Estimate_elem = SubElement(ConfigurationL2_Item, "Estimate" + chain_id)

        SpectralShiftFiltering = SubElement(Estimate_elem, "SpectralShiftFiltering")
        SpectralShiftFiltering.text = str(params_curr.spectral_shift_filtering)

        ProductResolution = SubElement(Estimate_elem, "ProductResolution")
        ProductResolution.text = str(params_curr.product_resolution)

        model_params = params_curr.model_parameters
        ModelParameters = SubElement(Estimate_elem, "ModelParameters")

        MaximumHeight = SubElement(ModelParameters, "MaximumHeight")
        MaximumHeight.text = str(model_params.maximum_height)

        NumberOfExtinctionValue = SubElement(ModelParameters, "NumberOfExtinctionValue")
        NumberOfExtinctionValue.text = str(model_params.number_of_extinction_value)

        NumberOfGroundVolumeRatioValue = SubElement(ModelParameters, "NumberOfGroundVolumeRatioValue")
        NumberOfGroundVolumeRatioValue.text = str(model_params.number_of_ground_volume_ratio_value)

        NumberOfTemporalDecorrelationValue = SubElement(ModelParameters, "NumberOfTemporalDecorrelationValue")
        NumberOfTemporalDecorrelationValue.text = str(model_params.number_of_temporal_decorrelation_value)

        MedianFactor = SubElement(Estimate_elem, "MedianFactor")
        MedianFactor.text = str(params_curr.median_factor)

    if chain_id == "FD":

        Estimate_elem = SubElement(ConfigurationL2_Item, "ChangeDetection")

        ProductResolution = SubElement(Estimate_elem, "ProductResolution")
        ProductResolution.text = str(params_curr.product_resolution)

        ConfidenceLevel = SubElement(Estimate_elem, "ConfidenceLevel")
        ConfidenceLevel.text = str(params_curr.confidence_level)

    if chain_id == "TOMO_FH" or chain_id == "TOMO":

        ProductResolution = SubElement(ConfigurationL2_Item, "ProductResolution")
        ProductResolution.text = str(params_curr.product_resolution)

        vertical_range = configuration_params.vertical_range

        VerticalRange = SubElement(ConfigurationL2_Item, "VerticalRange")

        MaximumHeight = SubElement(VerticalRange, "MaximumHeight")
        MaximumHeight.text = str(vertical_range.maximum_height)

        MinimumHeight = SubElement(VerticalRange, "MinimumHeight")
        MinimumHeight.text = str(vertical_range.minimum_height)

        sampling = SubElement(VerticalRange, "sampling")
        sampling.text = str(vertical_range.sampling)

    if chain_id == "TOMO_FH":

        Estimate_elem = SubElement(ConfigurationL2_Item, "EstimateFH")

        EnableSuperResolution = SubElement(Estimate_elem, "EnableSuperResolution")
        EnableSuperResolution.text = str(params_curr.enable_super_resolution)

        RegularizationNoiseFactor = SubElement(Estimate_elem, "RegularizationNoiseFactor")
        RegularizationNoiseFactor.text = str(params_curr.regularization_noise_factor)

        PowerThreshold = SubElement(Estimate_elem, "PowerThreshold")
        PowerThreshold.text = str(params_curr.power_threshold)

        MedianFactor = SubElement(Estimate_elem, "MedianFactor")
        MedianFactor.text = str(params_curr.median_factor)

    ConfigurationL2_Item = write_configuration_flags_core(ConfigurationL2_Item, configuration_params)

    # write to file
    ElementTree_indent(ConfigurationL2_Item)
    tree = ElementTree(ConfigurationL2_Item)
    tree.write(open(configuration_file_xml, "w"), encoding="unicode")


def write_agb_configuration_file(configuration_params, configuration_file_xml):
    """Write the configuration file of the AGB main APP"""

    chain_id = configuration_params.chain_id
    if not chain_id == "AGB":
        error_str = "This is not an AGB configuration file"
        logging.error(error_str)
        raise ValueError(error_str)

    ConfigurationL2_Item = Element("ConfigurationL2" + chain_id)
    ConfigurationL2_Item.set("version", _XML_VERSION)

    ground_canc_params = configuration_params.ground_cancellation
    ConfigurationL2_Item = write_ground_cancellation_core(ConfigurationL2_Item, ground_canc_params)

    Estimate_elem = SubElement(ConfigurationL2_Item, "EstimateAGB")
    write_estimateagb_core(Estimate_elem, configuration_params)

    ConfigurationL2_Item = write_configuration_flags_core(ConfigurationL2_Item, configuration_params)

    # write to file
    ElementTree_indent(ConfigurationL2_Item)
    tree = ElementTree(ConfigurationL2_Item)
    tree.write(open(configuration_file_xml, "w"), encoding="unicode")


def write_coreprocessing_agb_configuration_file(configuration_params, configuration_file_xml):
    """Write the configuration file of the Core Processing AGB APP"""

    ConfigurationL2_Item = Element("ConfigurationL2CoreProcessingAGB")
    ConfigurationL2_Item.set("version", _XML_VERSION)

    Estimate_elem = SubElement(ConfigurationL2_Item, "EstimateAGB")
    write_estimateagb_core(Estimate_elem, configuration_params)

    # write to file
    ElementTree_indent(ConfigurationL2_Item)
    tree = ElementTree(ConfigurationL2_Item)
    tree.write(open(configuration_file_xml, "w"), encoding="unicode")


# common functions used by the writers (core):


def write_configuration_flags_core(ConfigurationL2_Item, configuration_params):
    enable_resampling = SubElement(ConfigurationL2_Item, "EnableResampling")
    enable_resampling.text = str(configuration_params.enable_resampling)
    compute_geometry = SubElement(ConfigurationL2_Item, "ComputeGeometry")
    compute_geometry.text = str(configuration_params.compute_geometry)
    apply_calibration_screen = SubElement(ConfigurationL2_Item, "ApplyCalibrationScreen")
    apply_calibration_screen.text = str(configuration_params.apply_calibration_screen)
    DEM_flattening = SubElement(ConfigurationL2_Item, "DEMflattening")
    DEM_flattening.text = str(configuration_params.DEM_flattening)
    multilook_heading_correction = SubElement(ConfigurationL2_Item, "MultilookHeadingCorrection")
    multilook_heading_correction.text = str(configuration_params.multilook_heading_correction)
    save_breakpoints = SubElement(ConfigurationL2_Item, "SaveBreakpoints")
    save_breakpoints.text = str(configuration_params.save_breakpoints)
    delete_temporary_files = SubElement(ConfigurationL2_Item, "DeleteTemporaryFiles")
    delete_temporary_files.text = str(configuration_params.delete_temporary_files)

    return ConfigurationL2_Item


def write_ground_cancellation_core(ConfigurationL2_Item, ground_canc_params):

    GroundCancellation = SubElement(ConfigurationL2_Item, "GroundCancellation")
    MultiMaster = SubElement(GroundCancellation, "MultiMaster")
    MultiMaster.text = str(ground_canc_params.multi_master_flag)
    EnhancedForestHeight = SubElement(GroundCancellation, "EnhancedForestHeight")
    EnhancedForestHeight.text = str(ground_canc_params.enhanced_forest_height)
    EnhancedForestHeight.set("unit", "m")
    ModelBasedEqualization = SubElement(GroundCancellation, "ModelBasedEqualization")
    ModelBasedEqualization.text = str(ground_canc_params.equalization_flag)
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
    return ConfigurationL2_Item


def write_estimateagb_core(Estimate_elem, configuration_params):

    params_curr = getattr(configuration_params, "AGB")

    if params_curr.residual_function:
        Estimate_elem = write_residual_function_core(Estimate_elem, configuration_params)

    number_of_tests = SubElement(Estimate_elem, "number_of_tests")
    number_of_tests.text = str(params_curr.number_of_tests)

    forest_class_observable_name = SubElement(Estimate_elem, "forest_class_observable_name")
    forest_class_observable_name.text = str(params_curr.forest_class_observable_name)
    

    transfer_function_name = SubElement(Estimate_elem, "transfer_function_name")
    transfer_function_name.text = str(params_curr.transfer_function_name)
    
    fraction_of_roi_per_test = SubElement(Estimate_elem, "fraction_of_roi_per_test")
    fraction_of_roi_per_test.text = str(params_curr.fraction_of_roi_per_test)

    fraction_of_cal_per_test = SubElement(Estimate_elem, "fraction_of_cal_per_test")
    fraction_of_cal_per_test.text = str(params_curr.fraction_of_cal_per_test)

    add_variability_on_cal_data = SubElement(Estimate_elem, "add_variability_on_cal_data")
    add_variability_on_cal_data.text = str(params_curr.add_variability_on_cal_data)

    intermediate_ground_averaging = SubElement(Estimate_elem, "intermediate_ground_averaging")
    intermediate_ground_averaging.text = str(params_curr.intermediate_ground_averaging)

    product_resolution = SubElement(Estimate_elem, "product_resolution")
    product_resolution.text = str(params_curr.product_resolution)
    product_resolution.set("unit", "m")

    distance_sampling_area = SubElement(Estimate_elem, "distance_sampling_area")
    distance_sampling_area.text = str(params_curr.distance_sampling_area)
    distance_sampling_area.set("unit", "m")

    parameter_block_size = SubElement(Estimate_elem, "parameter_block_size")
    parameter_block_size.text = str(params_curr.parameter_block_size)
    parameter_block_size.set("unit", "m")

    distance_parameter_block = SubElement(Estimate_elem, "distance_parameter_block")
    distance_parameter_block.text = str(params_curr.distance_parameter_block)
    distance_parameter_block.set("unit", "m")

    min_number_of_rois = SubElement(Estimate_elem, "min_number_of_rois")
    min_number_of_rois.text = str(params_curr.min_number_of_rois)

    min_number_of_rois_per_stack = SubElement(Estimate_elem, "min_number_of_rois_per_stack")
    min_number_of_rois_per_stack.text = str(params_curr.min_number_of_rois_per_stack)

    min_number_of_cals_per_test = SubElement(Estimate_elem, "min_number_of_cals_per_test")
    min_number_of_cals_per_test.text = str(params_curr.min_number_of_cals_per_test)

    min_number_of_rois_per_test = SubElement(Estimate_elem, "min_number_of_rois_per_test")
    min_number_of_rois_per_test.text = str(params_curr.min_number_of_rois_per_test)

    estimation_valid_values_limits = SubElement(Estimate_elem, "EstimationValidValuesLimits")
    min_item = SubElement(estimation_valid_values_limits, "min")
    min_item.text = str(params_curr.estimation_valid_values_limits.min)
    min_item.set("unit", "t/ha")
    max_item = SubElement(estimation_valid_values_limits, "max")
    max_item.text = str(params_curr.estimation_valid_values_limits.max)
    max_item.set("unit", "t/ha")

    return Estimate_elem


def write_residual_function_core(Estimate_elem, configuration_params):

    residual_function = SubElement(Estimate_elem, "residual_function")
    
    formula_terms = SubElement(residual_function, "formula_terms")
    number_of_terms = len(configuration_params.AGB.residual_function.formula_terms.string)
    term_struct = configuration_params.AGB.residual_function.formula_terms
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
    number_of_parameters = len(configuration_params.AGB.residual_function.formula_parameters.name)
    par_struct = configuration_params.AGB.residual_function.formula_parameters
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

        associated_observable = SubElement(par,"associated_observable_name")
        associated_observable.text = par_struct.associated_observable_name[index]
        
        
        parameter_variabilities = SubElement(par,"variability")
        samples_item = SubElement(parameter_variabilities,"samples")
        samples_item.text = str(par_struct.parameter_variabilities[index][0])
        forest_classes_item = SubElement(parameter_variabilities,"forest_classes")
        forest_classes_item.text = str(par_struct.parameter_variabilities[index][1])
        stacks_item = SubElement(parameter_variabilities,"stacks")
        stacks_item.text = str(par_struct.parameter_variabilities[index][2])
        global_cycles_item = SubElement(parameter_variabilities,"global_cycles")
        global_cycles_item.text = str(par_struct.parameter_variabilities[index][3])
        headings_item = SubElement(parameter_variabilities,"headings")
        headings_item.text = str(par_struct.parameter_variabilities[index][4])
        swaths_item = SubElement(parameter_variabilities,"swaths")
        swaths_item.text = str(par_struct.parameter_variabilities[index][5])
        subswaths_item = SubElement(parameter_variabilities,"subswaths")
        subswaths_item.text = str(par_struct.parameter_variabilities[index][6])
        azimuth_images_item = SubElement(parameter_variabilities,"azimuth_images")
        azimuth_images_item.text = str(par_struct.parameter_variabilities[index][7])
        
        
    formula_observables = SubElement(residual_function, "formula_observables")
    obs_struct = configuration_params.AGB.residual_function.formula_observables
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
        sources.set("unit",str(obs_struct.source_unit[index_obs]))
        sources.set("resolution_m",str(obs_struct.source_resolution[index_obs]))
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

    return Estimate_elem


def write_fd_lookup_table(lookup_table_file_name_xml, table_type, lut_dict):
    """Write the lookup table for fd chain fnf and covariance  outputs"""

    # table_type: string 'fnf' or 'covariance'
    if not table_type == "fnf" and not table_type == "covariance":
        raise ValueError(
            '#3th iinput should be a string of value "fnf" or "covariance": "' + table_type + '" is not recognized'
        )

    if table_type == "fnf":
        lut_Item = Element("FDlookUpTableFNFmasks")
    elif table_type == "covariance":
        lut_Item = Element("FDlookUpTableCovariances")

    lut_Item.set("version", _XML_VERSION)

    for stack_id in lut_dict.keys():

        lookup_table_struct = lut_dict[stack_id]

        L1cStack = SubElement(lut_Item, "L1cStack")
        L1cStack.set("stack_id", stack_id)

        GeographicBoundaries = SubElement(L1cStack, "GeographicBoundaries")
        GeographicBoundaries.set("unit", "deg")

        lat_min = SubElement(GeographicBoundaries, "latMin")
        lat_min.text = lookup_table_struct.boundaries.lat_min
        lat_max = SubElement(GeographicBoundaries, "latMax")
        lat_max.text = lookup_table_struct.boundaries.lat_max
        lon_min = SubElement(GeographicBoundaries, "lonMin")
        lon_min.text = lookup_table_struct.boundaries.lon_min
        lon_max = SubElement(GeographicBoundaries, "lonMax")
        lon_max.text = lookup_table_struct.boundaries.lon_max

        for step_str_idx in lookup_table_struct.dates.keys():

            dates_and_fnames_element = SubElement(L1cStack, table_type)
            dates_and_fnames_element.set("step", step_str_idx)

            date = SubElement(dates_and_fnames_element, "Date")
            date.text = lookup_table_struct.dates[step_str_idx]
            date.set("unit", "Utc")

            file_name = SubElement(dates_and_fnames_element, "FileName")
            file_name.text = lookup_table_struct.file_names[step_str_idx]

    # write to file
    ElementTree_indent(lut_Item)
    tree = ElementTree(lut_Item)
    tree.write(open(lookup_table_file_name_xml, "w"), encoding="unicode")


def parse_biomassL2_main_input_file(input_file_xml):
    """Parse the input file, configuration file and current campaign params file
    and store all the needed information in the dataset data_stack class"""

    # 1) Parse input_file_xml -------------------------------------------------
    tree = ET.parse(input_file_xml)
    root = tree.getroot()

    L1c_repository = root.find("L1cRepository").text
    if not os.path.exists(L1c_repository):
        error_message = "Main input file: the specified L1c Repository folder does not exist "
        logging.error(error_message)
        raise ValueError(error_message)

    L1c_aux_data_repository = root.find("AuxiliaryProductsFolder").text
    if not os.path.exists(L1c_aux_data_repository):
        error_message = "Main input file: the specified AuxiliaryProductsFolder does not exist: "
        logging.error(error_message)
        raise ValueError(error_message)

    dem_folder_found = False
    aux_list = os.listdir(L1c_aux_data_repository)
    for aux_name in aux_list:
        if aux_name == "DEM":
            dem_folder_found = True
    if not dem_folder_found:
        error_message = (
            'Main input file: the specified AuxiliaryProductsFolder should contain AT LEAST the sub-folder "DEM"'
        )
        logging.error(error_message)
        raise ValueError(error_message)

    output_folder = root.find("OutputFolder").text

    try:
        geographic_grid_sampling = float(root.find("GeographicGridSampling").text)
    except:
        error_message = "Main input file: GeographicGridSampling should be a numeric value"
        logging.error(error_message)
        raise ValueError(error_message)
    if geographic_grid_sampling <= 0:
        error_message = "Main input file: GeographicGridSampling should be a positive numeric value [m]"
        logging.error(error_message)
        raise ValueError(error_message)

    L2Product_Item = root.find("L2Product")
    proc_flags_AGB = bool_from_string(L2Product_Item.find("AGB").text)
    proc_flags_FD = bool_from_string(L2Product_Item.find("FD").text)
    proc_flags_FH = bool_from_string(L2Product_Item.find("FH").text)
    proc_flags_TOMO = bool_from_string(L2Product_Item.find("TOMO").text)
    proc_flags_TOMO_FH = bool_from_string(L2Product_Item.find("TOMO_FH").text)

    proc_flags_struct = proc_flags(proc_flags_AGB, proc_flags_FH, proc_flags_TOMO_FH, proc_flags_FD, proc_flags_TOMO)

    L1cDates = root.findall("L1cDate")
    if len(L1cDates) != 2:
        error_message = 'Main input file: there should be #2 L1cDates (one with value="start" and one with value="stop"'
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

    point_list = []
    GeographicBoundary_Item = root.find("GeographicBoundary")
    Points = GeographicBoundary_Item.findall("Point")
    if len(Points) < 3:
        error_message = "Main input file: at least #3 Geograplic Boundary points should be specified"
        logging.error(error_message)
        raise ValueError(error_message)

    for Point_Item in Points:

        Latitude_str = Point_Item.find("Latitude").text
        Longitude_str = Point_Item.find("Longitude").text
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

    stacks_to_find_list = []
    Stacks_to_find_Item = root.find("Stacks_to_find")
    if Stacks_to_find_Item != None:
        stacks_to_find = Stacks_to_find_Item.findall("Stack_id")
        if len(stacks_to_find):
            for stack_curr in stacks_to_find:
                stacks_to_find_list.append(stack_curr.text)
                if stack_curr.text == None:
                    error_message = 'Main input file: "Stacks to find" optional input cannot have empty "Stack_id" fields: fill them or totally erase the "Stacks to find" element from input file'
                    logging.error(error_message)
                    raise ValueError(error_message)

    main_input_struct = main_input_params(
        L1c_repository,
        L1c_aux_data_repository,
        proc_flags_struct,
        L1c_date_struct,
        geographic_boundary_struct,
        output_folder,
        geographic_grid_sampling,
        stacks_to_find_list,
    )
    return main_input_struct


def parse_chains_input_file(input_file_xml):
    """Parse the input file of inner chains(AGB, FH, FD, TOMO FH and TOMO)"""

    # 1) Parse input_file_xml -------------------------------------------------
    tree = ET.parse(input_file_xml)
    root = tree.getroot()
    chain_id = root.tag[8:]
    if chain_id != "AGB" and chain_id != "FH" and chain_id != "FD" and chain_id != "TOMO_FH" and chain_id != "TOMO":
        error_message = "The provided input file is not an inner biomassL2 chain input (AGB, FH, FD, TOMO FH or TOMO)"
        logging.error(error_message)
        raise ValueError(error_message)

    L1cRepository_Item = root.find("L1cRepository")
    L1c_repository = L1cRepository_Item.text

    L1cProductList_Item = root.find("L1cProductList")
    L1cProduct_Items = L1cProductList_Item.findall("L1cProduct")
    number_of_stacks = len(L1cProduct_Items)
    stack_composition = {}
    for L1cProduct_Item in L1cProduct_Items:

        stack_id = L1cProduct_Item.attrib["unique_stack_id"]

        Acquisition_Items = L1cProduct_Item.findall("Acquisition")
        Acquisitions_list = []
        for Acquisition_Item in Acquisition_Items:
            Acquisitions_list.append(Acquisition_Item.text)

        stack_composition[stack_id] = Acquisitions_list

    AuxiliaryProductList_Item = root.find("AuxiliaryProductList")
    DEM_Item = AuxiliaryProductList_Item.find("DEM")
    dem_folder = DEM_Item.text

    if chain_id == "FD":
        AVG_Item = AuxiliaryProductList_Item.find("AverageCovarianceFolder")
        average_covariance_folder = AVG_Item.text
    else:
        average_covariance_folder = ""

    L1cStack_Items = AuxiliaryProductList_Item.findall("L1cStack")
    if len(L1cStack_Items) != number_of_stacks:
        error_message = (
            chain_id
            + " Input File: the number of L1cStack elements should be equal to the number of L1cProduct elements"
        )
        logging.error(error_message)
        raise ValueError(error_message)

    ECEF_grid_file_names = {}
    kz_file_names = {}
    off_nadir_angle_file_names = {}
    slope_file_names = {}
    reference_height_file_names = {}
    slant_range_distances_file_names = {}
    calibration_screens_file_names = {}
    for L1cStack_Item in L1cStack_Items:

        stack_id = L1cStack_Item.attrib["unique_stack_id"]
        if not stack_id in stack_composition.keys():
            error_message = (
                chain_id
                + ' Input File: the "stack_id" values of L1cStack elelemts does not match the "stack_id" values of L1cProduct elements'
            )
            logging.error(error_message)
            raise ValueError(error_message)

        Geometry_Item = L1cStack_Item.find("Geometry")

        ECEF_grid_file_names[stack_id] = Geometry_Item.find("ECEFGrid").text
        kz_file_names[stack_id] = Geometry_Item.find("kz").text
        slant_range_distances_file_names[stack_id] = Geometry_Item.find("RadarDistances").text
        off_nadir_angle_file_names[stack_id] = Geometry_Item.find("OffNadirAngles").text
        slope_file_names[stack_id] = Geometry_Item.find("Slope").text
        reference_height_file_names[stack_id] = Geometry_Item.find("ReferenceHeight").text

        if L1cStack_Item.find("CalibrationScreens") != None:
            calibration_screens_file_names[stack_id] = L1cStack_Item.find("CalibrationScreens").text
        elif chain_id == "TOMO" or chain_id == "TOMO_FH":
            error_message = chain_id + " Input File: CalibrationScreens are mandatory for " + chain_id + " chain."
            logging.error(error_message)
            raise ValueError(error_message)

    reference_agb_folder = ""
    system_decorrelation_fun_folder = ""
    forest_mask_catalogue_folder = ""
    forest_height_folder = ""

    if chain_id == "AGB" and AuxiliaryProductList_Item.find("ReferenceAGB").text:
        reference_agb_folder = AuxiliaryProductList_Item.find("ReferenceAGB").text
    if chain_id == "AGB" and AuxiliaryProductList_Item.find("ForestHeight").text:
        forest_height_folder = AuxiliaryProductList_Item.find("ForestHeight").text
    if chain_id == "FH":
        if AuxiliaryProductList_Item.find("SystemDecorrelationFunction"):
            system_decorrelation_fun_folder = AuxiliaryProductList_Item.find("SystemDecorrelationFunction").text
        else:
            system_decorrelation_fun_folder = ""

    if not "TOMO" in chain_id:
        forest_mask_catalogue_folder = AuxiliaryProductList_Item.find("ForestMask").text

    output_folder = root.find("OutputFolder").text

    geographic_grid_sampling = float(root.find("GeographicGridSampling").text)

    proc_inputs = input_params(
        chain_id,
        L1c_repository,
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
        output_folder,
        geographic_grid_sampling,
        system_decorrelation_fun_folder,
        average_covariance_folder,
        forest_height_folder,
    )

    return proc_inputs


def parse_fd_lookup_table(lookup_table_file_name_xml):
    """Parse the look up tables uded to indicize the fnf mask and covariance at each FD step"""

    lut_dict = {}
    # 1) Parse input_file_xml -------------------------------------------------
    tree = ET.parse(lookup_table_file_name_xml)
    root = tree.getroot()

    L1cStacks = root.findall("L1cStack")
    for L1cStack in L1cStacks:

        lookup_table_struct = namedtuple("lookup_table", "boundaries dates file_names")
        stack_id = L1cStack.attrib["stack_id"]

        # extract all the bondaries coordinates from current stack (one date for stack):
        geographic_boundaries_curr = geographic_boundaries(None, None, None, None)

        GeographicBoundaries_Item = L1cStack.find("GeographicBoundaries")
        geographic_boundaries_curr.lat_min = GeographicBoundaries_Item.find("latMin").text
        geographic_boundaries_curr.lat_max = GeographicBoundaries_Item.find("latMax").text
        geographic_boundaries_curr.lon_min = GeographicBoundaries_Item.find("lonMin").text
        geographic_boundaries_curr.lon_max = GeographicBoundaries_Item.find("lonMax").text

        # extract all the dates from current stack (one date for each step):
        dates_dict = {}
        fnames_dict = {}
        try:
            if L1cStack.getchildren()[0].tag == "GeographicBoundaries":
                table_type = L1cStack.getchildren()[1].tag
            else:
                table_type = L1cStack.getchildren()[0].tag
        except:
            raise ImportError(" table xml is not valid")

        if not table_type == "fnf" and not table_type == "covariance":
            raise ImportError(' table xml is not valid, element "' + table_type + '" not recognized')

        fd_steps = L1cStack.findall(table_type)
        for step_curr in fd_steps:
            step_number = step_curr.attrib["step"]
            dates_dict[step_number] = step_curr.find("Date").text
            fnames_dict[step_number] = step_curr.find("FileName").text

        # fill the optupt dictionary with the current stack boundaries and dates:
        lookup_table_struct.boundaries = geographic_boundaries_curr
        lookup_table_struct.dates = dates_dict
        lookup_table_struct.file_names = fnames_dict
        lut_dict[stack_id] = lookup_table_struct

        del lookup_table_struct, dates_dict

    return lut_dict


def parse_biopal_configuration_file(biopal_configuration_file_xml):
    """Parse the biopal configuration file"""

    gdal_path = gdal_environment_path = None

    tree = ET.parse(biopal_configuration_file_xml)
    root = tree.getroot()

    gdal_item = root.find("GDAL")
    if gdal_item:
        gdal_path = gdal_item.find("gdal_path").text
        gdal_environment_path = gdal_item.find("gdal_environment_path").text

    return gdal_path, gdal_environment_path


def parse_chains_configuration_file(configuration_file_xml, output_folder=""):
    """Parse the chains configuration files for FD, FH, TOMO, TOMO_FH, not AGB"""

    tree = ET.parse(configuration_file_xml)
    root = tree.getroot()
    chain_id = root.tag[15:]

    if chain_id == "AGB":
        raise ValueError("Please, use the parse_agb_configuration_file and parse_coreprocessing_agb_configuration_file")

    # interpolate_stack
    if chain_id == "TOMO":
        InterpolateStack_Item = root.find("InterpolateStack")
        regular_baseline_grid_flag = bool_from_string(InterpolateStack_Item.find("RegularBaselineGrid").text)
        interpolate_stack = interpolate_params(regular_baseline_grid_flag)
    else:
        interpolate_stack = None

    # ground_cancellation
    ground_cancellation = None
    if chain_id == "FD":
        ground_cancellation = parse_ground_cancellation_core(root)

    # vertical_range
    if chain_id == "TOMO" or chain_id == "TOMO_FH":
        VerticalRange_Item = root.find("VerticalRange")
        maximum_height = float(VerticalRange_Item.find("MaximumHeight").text)
        minimum_height = float(VerticalRange_Item.find("MinimumHeight").text)
        sampling = float(VerticalRange_Item.find("sampling").text)
        vertical_range = vertical_range_params(maximum_height, minimum_height, sampling)
    else:
        vertical_range = None

    (
        enable_resampling,
        compute_geometry,
        apply_calibration_screen,
        DEM_flattening,
        multilook_heading_correction,
        save_breakpoints,
        delete_temporary_files,
    ) = parse_configuration_flags_core(root)

    if chain_id == "FD":
        chain_field_name = "ChangeDetection"
    elif chain_id == "FH" or chain_id == "TOMO_FH":
        chain_field_name = "EstimateFH"
    elif chain_id == "TOMO":
        chain_field_name = ""

    chain_field_Item = root.find(chain_field_name)

    if chain_id == "TOMO" or chain_id == "TOMO_FH":
        product_resolution = float(root.find("ProductResolution").text)
    else:
        product_resolution = float(chain_field_Item.find("ProductResolution").text)

    AGB = None
    FH = None
    FD = None
    TOMO_FH = None
    TOMO = None
    if chain_id == "FD":

        confidence_level = float(chain_field_Item.find("ConfidenceLevel").text)
        FD = FD_est_params(product_resolution, confidence_level)

    elif chain_id == "FH":

        estimation_valid_values_limits = [
            float(chain_field_Item.find("EstimationValidValuesLimits").find("min").text),
            float(chain_field_Item.find("EstimationValidValuesLimits").find("max").text),
        ]

        spectral_shift_filtering = bool_from_string(chain_field_Item.find("SpectralShiftFiltering").text)

        kz_thresholds = [
            float(chain_field_Item.find("KZ_thresholds").find("min").text),
            float(chain_field_Item.find("KZ_thresholds").find("max").text),
        ]

        median_factor_temp = float(chain_field_Item.find("MedianFactor").text)
        median_factor = int(median_factor_temp)
        if median_factor_temp - median_factor > 0:
            logging.warning(
                "FH configuration file MedianFactor should be integer, the value that will be used is: {}".format(
                    median_factor
                )
            )

        model_parameters_Item = chain_field_Item.find("ModelParameters")
        maximum_height = int(model_parameters_Item.find("MaximumHeight").text)
        number_of_extinction_value = int(model_parameters_Item.find("NumberOfExtinctionValue").text)
        number_of_ground_volume_ratio_value = int(model_parameters_Item.find("NumberOfGroundVolumeRatioValue").text)
        number_of_temporal_decorrelation_value = int(
            model_parameters_Item.find("NumberOfTemporalDecorrelationValue").text
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
        model_parameters = FH_model_params(
            estimation_valid_values_limits,
            maximum_height,
            number_of_extinction_value,
            number_of_ground_volume_ratio_value,
            number_of_temporal_decorrelation_value,
        )

        FH = FH_est_params(
            spectral_shift_filtering, product_resolution, kz_thresholds, model_parameters, median_factor,
        )

    elif chain_id == "TOMO_FH":

        estimation_valid_values_limits = [
            float(chain_field_Item.find("EstimationValidValuesLimits").find("min").text),
            float(chain_field_Item.find("EstimationValidValuesLimits").find("max").text),
        ]

        enable_super_resolution = bool_from_string(chain_field_Item.find("EnableSuperResolution").text)
        regularization_noise_factor = float(chain_field_Item.find("RegularizationNoiseFactor").text)
        power_threshold = float(chain_field_Item.find("PowerThreshold").text)
        median_factor = int(chain_field_Item.find("MedianFactor").text)
        TOMO_FH = TOMO_FH_est_params(
            estimation_valid_values_limits,
            product_resolution,
            enable_super_resolution,
            regularization_noise_factor,
            power_threshold,
            median_factor,
        )

    elif chain_id == "TOMO":

        TOMO = TOMO_est_params(product_resolution)

    proc_config = configuration_params(
        chain_id,
        interpolate_stack,
        ground_cancellation,
        vertical_range,
        AGB,
        FH,
        FD,
        TOMO_FH,
        TOMO,
        enable_resampling,
        compute_geometry,
        apply_calibration_screen,
        DEM_flattening,
        multilook_heading_correction,
        save_breakpoints,
        delete_temporary_files,
    )

    return proc_config


def parse_agb_configuration_file(configuration_file_xml):
    "Parse the chains configuration files for AGB Processor (main APP)"

    tree = ET.parse(configuration_file_xml)
    root = tree.getroot()
    chain_id = root.tag[15:]

    if not chain_id == "AGB":
        error_str = "This is not an AGB configuration file"
        logging.error(error_str)

    interpolate_stack = None
    vertical_range = None

    ground_cancellation = parse_ground_cancellation_core(root)

    chain_field_Item = root.find("EstimateAGB")

    AGB = parse_estimateagb_core(chain_field_Item)

    (
        enable_resampling,
        compute_geometry,
        apply_calibration_screen,
        DEM_flattening,
        multilook_heading_correction,
        save_breakpoints,
        delete_temporary_files,
    ) = parse_configuration_flags_core(root)

    proc_config = configuration_params(
        chain_id,
        interpolate_stack,
        ground_cancellation,
        vertical_range,
        AGB,
        None,
        None,
        None,
        None,
        enable_resampling,
        compute_geometry,
        apply_calibration_screen,
        DEM_flattening,
        multilook_heading_correction,
        save_breakpoints,
        delete_temporary_files,
    )

    return proc_config


def parse_coreprocessing_agb_configuration_file(configuration_file_xml, output_folder=""):
    "Parse the configuration files for CoreProcessingAGB APP"

    tree = ET.parse(configuration_file_xml)
    root = tree.getroot()
    chain_id = root.tag[15:]

    # interpolate_stack
    if not chain_id == "CoreProcessingAGB":
        error_str = "This is not a Core Processing AGB APP configuration file"
        logging.error(error_str)
        raise ValueError(error_str)

    interpolate_stack = None
    vertical_range = None
    enable_resampling = None
    compute_geometry = None
    apply_calibration_screen = None
    DEM_flattening = None
    multilook_heading_correction = None
    save_breakpoints = None
    delete_temporary_files = None

    ground_cancellation = parse_ground_cancellation_core(root)

    chain_field_Item = root.find("EstimateAGB")

    AGB = parse_estimateagb_core(chain_field_Item, output_folder)

    chain_id_final = "AGB"

    proc_config = configuration_params(
        chain_id_final,
        interpolate_stack,
        ground_cancellation,
        vertical_range,
        AGB,
        None,
        None,
        None,
        None,
        enable_resampling,
        compute_geometry,
        apply_calibration_screen,
        DEM_flattening,
        multilook_heading_correction,
        save_breakpoints,
        delete_temporary_files,
    )

    return proc_config


def parse_ground_cancellation_core(root):
    ground_cancellation = None
    GroundCancellation_Item = root.find("GroundCancellation")
    if GroundCancellation_Item:
        multi_master_flag = bool_from_string(GroundCancellation_Item.find("MultiMaster").text)
        enhanced_forest_height = float(GroundCancellation_Item.find("EnhancedForestHeight").text)
        equalization_flag = GroundCancellation_Item.find("ModelBasedEqualization").text
        if not equalization_flag == "1" and not equalization_flag == "2" and not equalization_flag == "3":
            error_str = 'Configuration flag "ModelBasedEqualization" value "{}" not valid, choose among "1", "2" or "3", where "1"->always OFF; "2"->always ON; "3"->ON only if two acquisitions are present'.format(
                equalization_flag
            )
            logging.error(error_str)
            raise ValueError(error_str)
        ground_cancellation = ground_canc_params(multi_master_flag, enhanced_forest_height, equalization_flag)
    return ground_cancellation


def parse_configuration_flags_core(root):
    enable_resampling = bool_from_string(root.find("EnableResampling").text)
    compute_geometry = bool_from_string(root.find("ComputeGeometry").text)
    apply_calibration_screen = bool_from_string(root.find("ApplyCalibrationScreen").text)
    DEM_flattening = bool_from_string(root.find("DEMflattening").text)
    multilook_heading_correction = bool_from_string(root.find("MultilookHeadingCorrection").text)
    save_breakpoints = bool_from_string(root.find("SaveBreakpoints").text)
    delete_temporary_files = bool_from_string(root.find("DeleteTemporaryFiles").text)

    return (
        enable_resampling,
        compute_geometry,
        apply_calibration_screen,
        DEM_flattening,
        multilook_heading_correction,
        save_breakpoints,
        delete_temporary_files,
    )


def parse_estimateagb_core(chain_field_Item, output_folder=""):
    residual_function_struct = None
    residual_function_Item = chain_field_Item.find("residual_function")
    if residual_function_Item:
        residual_function_struct = parse_agb_residual_function_core(residual_function_Item, output_folder)

    number_of_tests = int(chain_field_Item.find("number_of_tests").text)
    forest_class_observable_name = str(chain_field_Item.find("forest_class_observable_name").text)
    transfer_function_name = str(chain_field_Item.find("transfer_function_name").text)
    intermediate_ground_averaging = float(chain_field_Item.find("intermediate_ground_averaging").text)

    
    distance_sampling_area = float(chain_field_Item.find("distance_sampling_area").text)
    fraction_of_roi_per_test = float(chain_field_Item.find("fraction_of_roi_per_test").text)
    fraction_of_cal_per_test = float(chain_field_Item.find("fraction_of_cal_per_test").text)

    add_variability_on_cal_data = bool_from_string(chain_field_Item.find("add_variability_on_cal_data").text)

    parameter_block_size = float(chain_field_Item.find("parameter_block_size").text)
    distance_parameter_block = float(chain_field_Item.find("distance_parameter_block").text)
    min_number_of_rois = int(chain_field_Item.find("min_number_of_rois").text)
    min_number_of_rois_per_stack = int(chain_field_Item.find("min_number_of_rois_per_stack").text)
    min_number_of_cals_per_test = int(chain_field_Item.find("min_number_of_cals_per_test").text)
    min_number_of_rois_per_test = int(chain_field_Item.find("min_number_of_rois_per_test").text)

    estimation_min_value = float(chain_field_Item.find("EstimationValidValuesLimits").find("min").text)
    estimation_max_value = float(chain_field_Item.find("EstimationValidValuesLimits").find("max").text)
    estimation_valid_values_limits = min_max(estimation_min_value, estimation_max_value,)

    product_resolution = float(chain_field_Item.find("product_resolution").text)

    AGB = agb_est_params(
        residual_function_struct,
        number_of_tests,
        forest_class_observable_name,
        transfer_function_name,
        fraction_of_roi_per_test,
        fraction_of_cal_per_test,
        add_variability_on_cal_data,
        intermediate_ground_averaging,
        product_resolution,
        distance_sampling_area,
        parameter_block_size,
        distance_parameter_block,
        min_number_of_rois,
        min_number_of_rois_per_stack,
        min_number_of_cals_per_test,
        min_number_of_rois_per_test,
        estimation_valid_values_limits,
    )
    return AGB


def parse_agb_residual_function_core(residual_function_Item, output_folder=""):

    formula_Item = residual_function_Item.find("formula_terms")
    formula_parameters_Item = residual_function_Item.find("formula_parameters")
    formula_observables_Item = residual_function_Item.find("formula_observables")

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
    formula_weight_struct = formula_weights(formula_weight_step1,formula_weight_step2,formula_weight_step3)
    formula_terms_struct = formula_terms(
        formula_name,
        formula_string,
        formula_weight_struct,
    )


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
        variabilities.append([
            bool_from_string(par_item.find("variability").find("samples").text),
            bool_from_string(par_item.find("variability").find("forest_classes").text),
            bool_from_string(par_item.find("variability").find("stacks").text),
            bool_from_string(par_item.find("variability").find("global_cycles").text),
            bool_from_string(par_item.find("variability").find("headings").text),
            bool_from_string(par_item.find("variability").find("swaths").text),
            bool_from_string(par_item.find("variability").find("subswaths").text),
            bool_from_string(par_item.find("variability").find("azimuth_images").text),
        ])
        
        
    formula_parameters_struct = formula_parameters(
        name,
        save_as_map,
        transform,
        limits,
        units_par,
        associated_observable,
        variabilities,
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

                    layer_path = os.path.join(os.path.dirname(output_folder), layer_item.text)
                    band_id = int(layer_item.attrib["band"])
                    curr_file_layers.append(layer_path)
                    curr_file_layers.append(band_id)
                    
                curr_stack_files.append(curr_file_layers)
            curr_obs_stacks.append(curr_stack_files)
        source_paths.append(curr_obs_stacks)

    formula_observables_struct = formula_observables(
        name, is_required, source_paths, source_units, source_resolutions, limits_units, limits, transform, averaging_method,
    )

    residual_function_struct = residual_function(formula_terms_struct, formula_parameters_struct, formula_observables_struct,)
    return residual_function_struct


def check_chains_input_file(proc_inputs):

    if not os.path.exists(proc_inputs.L1c_repository):
        error_message = " Input file L1cRepository does not exist: " + proc_inputs.L1c_repository
        logging.error(error_message)
        raise RuntimeError(error_message)

    for key in proc_inputs.stack_composition.keys():
        for pf_name in proc_inputs.stack_composition[key]:
            fullPath = os.path.join(proc_inputs.L1c_repository, pf_name)
            if not os.path.exists(fullPath):
                error_message = " Input file Acquisition " + pf_name + " does not exist: " + fullPath
                logging.error(error_message)
                raise RuntimeError(error_message)

    if not os.path.exists(proc_inputs.dem_folder):
        error_message = " Input file AuxiliaryProductList DEM does not exist: " + proc_inputs.dem_folder
        logging.error(error_message)
        raise RuntimeError(error_message)

    for ecef_name in proc_inputs.ECEF_grid_file_names:
        if not os.path.exists(ecef_name):
            error_message = " Input file AuxiliaryProductList ECEFgrid does not exist: " + ecef_name
            logging.error(error_message)
            raise RuntimeError(error_message)

    for kz_name in proc_inputs.kz_file_names:
        if not os.path.exists(kz_name):
            error_message = " Input file AuxiliaryProductList KZ does not exist: " + kz_name
            logging.error(error_message)
            raise RuntimeError(error_message)

    for inc_name in proc_inputs.off_nadir_angle_file_names:
        if not os.path.exists(inc_name):
            error_message = " Input file AuxiliaryProductList OffNadirAngle does not exist: " + inc_name
            logging.error(error_message)
            raise RuntimeError(error_message)

    for h_name in proc_inputs.reference_height_file_names:
        if not os.path.exists(h_name):
            error_message = "Input file AuxiliaryProductList ReferenceHeight does not exist: " + h_name
            logging.error(error_message)
            raise RuntimeError(error_message)

    chain_id = proc_inputs.chain_id
    if chain_id == "AGB":
        if not os.path.exists(proc_inputs.reference_agb_folder):
            error_message = (
                " Input file AuxiliaryProductList ReferenceAGB does not exist: " + proc_inputs.reference_agb_folder
            )
            logging.error(error_message)
            raise RuntimeError(error_message)
    if chain_id == "FD":
        if not os.path.exists(proc_inputs.average_covariance_folder):
            error_message = (
                " Input  AuxiliaryProductList AverageCovarianceFolder does not exist: "
                + proc_inputs.average_covariance_folder
            )
            logging.error(error_message)
            raise RuntimeError(error_message)
    if chain_id == "FH":
        if not os.path.exists(proc_inputs.system_decorrelation_fun_folder):
            error_message = (
                " Input file AuxiliaryProductList SystemDecorrelationFunction does not exist: "
                + proc_inputs.system_decorrelation_fun_folder
            )
            logging.error(error_message)
            raise RuntimeError(error_message)

    if not os.path.exists(proc_inputs.forest_mask_catalogue_folder):
        error_message = (
            " Input file AuxiliaryProductList forest non forest mask catalogue folder does not exist: "
            + proc_inputs.forest_mask_catalogue_folder
        )
        logging.error(error_message)
        raise RuntimeError(error_message)

    if len(proc_inputs.geographic_grid_sampling) == 0 or proc_inputs.geographic_grid_sampling <= 0:
        error_message = (
            " Input file GeographicGridSampling should be a positive number: " + proc_inputs.geographic_grid_sampling
        )
        logging.error(error_message)
        raise RuntimeError(error_message)

    if not os.path.exists(proc_inputs.forest_mask_catalogue_folder):
        error_message = (
            " Input file AuxiliaryProductList ForestMask does not exist: " + proc_inputs.forest_mask_catalogue_folder
        )
        logging.error(error_message)
        raise RuntimeError(error_message)


class XmlIO:
    """
    Class to read/write XML documents
    Allows python object-like access to parameters.
    pp = XmlIO('path_to_pp.xml')
    print(pp.v0)     # v0 is automatically converted to floating point
    print(pp.r[:10]) # first 10 values of range vector, also floating point
    The parameter object also allows dictionary-like access to handle problematic parameter names
    (which clash with python keywords). For example:
    print(pp['lambda']) # pp.lambda would be a syntax error
    print(pp['pass'])   # same as above
    print(pp['v0'])     # dictionary style access works for other parameters, too!
    The class provides full read/write support. Parameter values are changed by standard assignment
    and structures can be saved using the write method:
    pp.v0 = 100
    pp.write('path_to_new_pp.xml')
    """

    def __init__(self, root):
        if isinstance(root, str):
            self.__dict__["__root__"] = xml_tools.parse(root).find("object")
        else:
            self.__dict__["__root__"] = root

        if self.__root__ is None:
            raise ValueError('Expected an "object" element below the root element!')

        self.__dict__["__iterend__"] = False

    def __getstate__(self):
        return self.__root__

    def __setstate__(self, root):
        self.__dict__["__root__"] = root
        self.__dict__["__iterend__"] = False

    def __getparam__(self, name):
        p = [p for p in self.__root__.iter("parameter") if p.attrib["name"] == name]
        if len(p) != 1:
            raise AttributeError('Expected a unique match parameter name "%s", got %i matches.' % (name, len(p)))

        return [p[0].find(tag) for tag in ("remark", "datatype", "value", "unit")]

    @staticmethod
    def xml2val(v, t):
        type = t.text
        shape = t.attrib["length"]
        shape = np.asarray([np.uint64(d) for d in shape.split()])[::-1]
        size = np.prod(shape)

        if type == "pointer":
            p = v.find("parameter")
            return XmlIO.xml2val(*[p.find(t) for t in ("value", "datatype")])

        if type == "struct":
            obj_arr = [XmlIO(obj) for obj in v.iter("object")]
            return obj_arr[0] if size <= 1 else obj_arr

        conv = {
            "bool": bool,
            "int": int,
            "long": int,
            "float": np.float,
            "double": np.double,
            "string": lambda s: s,
        }
        try:
            if size > 1:
                val = np.asarray([conv[type](v) for v in v.text.strip("[]").split(",")]).reshape(shape)
            else:
                val = conv[type](v.text)
        except KeyError:
            print('XmlIO WARNING: Unsupported data type "%s" encountered. Skipping!' % (type))
            return None

        return val

    @staticmethod
    def val2xml(v, t, value):
        cdict = {
            str: (str, "string"),
            int: (str, "long"),
            float: (str, "double"),
            complex: (lambda z: "({},{})".format(z.real, z.imag), "complex"),
        }
        cdict[np.uint8] = cdict[int]
        cdict[np.int32] = cdict[int]
        cdict[np.uint32] = cdict[int]
        cdict[np.int64] = cdict[int]
        cdict[np.uint64] = cdict[int]
        cdict[np.float32] = (str, "float")
        cdict[np.float64] = cdict[float]
        cdict[np.complex64] = cdict[complex]
        cdict[bool] = cdict[str]

        if t.text == "pointer":
            p = v.find("parameter")
            return XmlIO.val2xml(*([p.find(t) for t in ("value", "datatype")] + [value]))

        try:
            vsize = 1 if isinstance(value, str) else len(value)
        except TypeError:
            vsize = 1

        t.attrib["length"] = str(vsize)
        v.clear()
        if vsize == 1 and not isinstance(value, XmlIO):
            t.text = cdict[type(value)][1]
            v.text = cdict[type(value)][0](value)
        elif all([isinstance(v, XmlIO) for v in value]):
            t.text = "struct"
            for obj in value:
                v.append(copy.deepcopy(obj.__root__))
        else:
            if isinstance(value, np.ndarray):
                t.attrib["length"] = " ".join([str(l) for l in value.shape[::-1]])
                value = value.flat
            vtype = type(value[0])
            t.text = cdict[vtype][1]
            v.text = "[" + ", ".join([cdict[vtype][0](val) for val in value]) + "]"

    def __getattr__(self, key):
        if key in self.__dict__:
            return self.__dict__[key]
        if key == 0:
            return self
        r, t, v, u = self.__getparam__(key)
        return XmlIO.xml2val(v, t)

    def __getitem__(self, key):
        return self.__getattr__(key)

    def __setattr__(self, name, value):
        if name in self.__dict__:
            self.__dict__[name] = value
            return
        r, t, v, u = self.__getparam__(name)
        XmlIO.val2xml(v, t, value)

    def __setitem__(self, key, value):
        self.__setattr__(key, value)

    def __contains__(self, key):
        try:
            _ = self.__getparam__(key)
        except AttributeError:
            return False
        return True

    def __len__(self):
        return 1

    def __iter__(self):
        self.__iterend__ = False
        return self

    def __next__(self):
        if self.__iterend__:
            raise StopIteration()
        self.__iterend__ = True
        return self

    def update(self, obj):
        try:
            d = obj.__dict__
        except AttributeError:
            d = obj
        if not isinstance(d, dict):
            raise ValueError("Expected a dictionary or an object with a __dict__ attribute!")

        for k in d:
            self.__setattr__(k, d[k])

    def __totree(self):
        ste_root = xml_tools.Element("stexml")
        ste_root.text = "\n"
        ste_root.append(copy.deepcopy(self.__root__))
        ste_root.addprevious(xml_tools.PI("xml-stylesheet", 'type="text/xsl" href="stexml.xsl"'))
        tree = xml_tools.ElementTree(ste_root)
        return tree

    def write(self, filename):
        self.__totree().write(filename, pretty_print=True, encoding="UTF-8", xml_declaration=True)

    def tostring(self):
        return xml_tools.tostring(self.__totree().getroot(), encoding="UTF-8")


def add_param(
    root, name, unit_text, datatype_text, remark_text="none", value_text=None, length=1, sub_flag=False,
):

    # Sub flag is needed in order to correctly add parameters to a struct-type xml tree
    if sub_flag:
        # Check if the value and object subelements were created previously
        if root.find("value") is None:
            value = ET.SubElement(root, "value")
            object = ET.SubElement(value, "object")
        else:
            root = root.find("value")
            object = root.find("object")
        root = object

    param = ET.SubElement(root, "parameter")
    param.set("name", name)

    unit = ET.SubElement(param, "unit")
    unit.text = unit_text

    datatype = ET.SubElement(param, "datatype")
    datatype.set("length", str(length))
    datatype.text = datatype_text

    remark = ET.SubElement(param, "remark")
    remark.text = remark_text

    if datatype_text != "struct":
        value = ET.SubElement(param, "value")
        value.text = value_text

    return param
