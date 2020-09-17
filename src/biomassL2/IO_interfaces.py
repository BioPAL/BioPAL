import struct
import os
import xml.etree.ElementTree as ET
from xml.etree.ElementTree import ElementTree, Element, SubElement
from collections import namedtuple
from namedlist import namedlist
from arepytools.timing.precisedatetime import PreciseDateTime
import ast
import logging

###############################################################################
# Fields to write and read the Biomass Header in the Binary Files:#############

_STRING_ENCODING = 'utf-8'
_date_str_STR_LEN = 33
_date_str_FORMAT = '{}s'.format(_date_str_STR_LEN)
_GEO_CORNERS_FORMAT = 'ffff'  # lon_min, lon_max, lat_min, lat_max
_UNIQUE_ACQ_ID_STR_LEN = 47
_UNIQUE_ACQ_ID_FORMAT = '{}s'.format(_UNIQUE_ACQ_ID_STR_LEN)
_HEADING_FORMAT = 'f'  # heading [deg]T
_HEADER_FORMAT = _date_str_FORMAT + _GEO_CORNERS_FORMAT + _UNIQUE_ACQ_ID_FORMAT
_HEADER_FORMAT_SIZE = struct.calcsize(_HEADER_FORMAT)
###############################################################################


###############################################################################
# structures (namedtuples) for read and write inputs and configuration files ##
_XML_VERSION = '1.0'


geographic_boundaries = namedlist(
    'geographic_boundaries',
    'lon_min \
                     lon_max \
                     lat_min \
                     lat_max',
)

# initialize output structure fields
raster_info = namedtuple(
    'raster_info',
    'num_samples \
 num_lines \
 pixel_spacing_slant_rg \
 pixel_spacing_az \
 resolution_m_slant_rg \
 resolution_m_az \
 lines_start_utc',
)

# Main Input parameters structure (namedtuple)
main_input_params = namedtuple(
    'main_input_params',
    'L1c_repository \
 L1c_aux_data_repository \
 proc_flags \
 L1c_date \
 geographic_boundary \
 output_folder \
 geographic_grid_sampling \
 stacks_to_find',
)

# main_input_params "proc_flags" sub-fields:
proc_flags = namedlist(
    'proc_flags_params',
    'AGB \
                         FH \
                         TOMO_FH \
                         FD \
                         TOMO',
)

# main_input_params "L1c_date" sub-fields:
L1c_date = namedtuple(
    'L1c_date_params',
    'start \
                       stop',
)

# main_input_params "L1c_date" sub-fields:
geographic_boundary = namedtuple('geographic_boundary_params', 'point')

# Input parameters structure (namedtuple)
input_params = namedtuple(
    'input_params',
    'chain_id \
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
 forest_height_folder',
)

# Configuration parameters structure (namedtuple)
configuration_params = namedtuple(
    'configuration_params',
    'chain_id \
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
 delete_temporary_files',
)

# define each sub-structure of the configuration parameters, with sub-namedtuples
# configuration_params "interpolate_stack"  sub-fields:
interpolate_params = namedtuple('interpolate_params', 'regular_baseline_grid_flag')
# configuration_params "ground_cancellation"  sub-fields:
ground_canc_params = namedtuple('ground_canc_params', 'multi_master_flag enhanced_forest_height equalization_flag')

# configuration_params "agb" sub-fields:
agb_est_params = namedtuple(
    'agb_est_params',
    'number_of_tests \
                                 fraction_of_roi_per_test \
                                 fraction_of_cal_per_test \
                                 intermediate_ground_averaging \
                                 product_resolution \
                                 distance_sampling_area \
                                 parameter_block_size \
                                 distance_parameter_block \
                                 min_number_of_rois \
                                 min_number_of_rois_per_stack \
                                 min_number_of_cals_per_test \
                                 min_number_of_rois_per_test \
                                 model_parameters \
                                 estimation_valid_values_limits',
)
# agb "model_parameters" sub-fields:
agb_model_params = namedtuple(
    'agb_model_params',
    'agb_model_scaling_limits agb_model_exponent_limits agb_cosine_exponent_limits agb_estimation_limits \
                              changes_across_pol changes_across_stacks changes_across_global_cycle changes_across_heading \
                              changes_across_swath changes_across_subswath changes_across_azimuth',
)
# model_parameters "triplet_params" sub-fields:
triplet_params = namedtuple('triplet_params', 'agb_model_scaling agb_model_exponent agb_cosine_exponent')

# configuration_params "FH" sub-fields:
FH_est_params = namedtuple(
    'FH_est_params',
    'spectral_shift_filtering \
                                product_resolution \
                                kz_thresholds \
                                model_parameters \
                                median_factor',
)
# FH "model_parameters" sub-fields:
FH_model_params = namedtuple(
    'FH_model_params',
    'estimation_valid_values_limits maximum_height number_of_extinction_value number_of_ground_volume_ratio_value number_of_temporal_decorrelation_value',
)

# configuration_params "FD" sub-fields:
FD_est_params = namedtuple(
    'FD_est_params',
    'product_resolution \
                                confidence_level',
)

# configuration_params "TOMO FH" sub-fields:
TOMO_FH_est_params = namedtuple(
    'TOMO_FH_est_params',
    'estimation_valid_values_limits \
                                 product_resolution \
                                 enable_super_resolution \
                                 regularization_noise_factor \
                                 power_threshold \
                                 median_factor',
)

# configuration_params "TOMO" sub-fields:
TOMO_est_params = namedtuple('TOMO_est_params', 'product_resolution')

# TOMO and TOMO_FH "vertical range" sub-fields:
vertical_range_params = namedtuple('vertical_range_params', 'maximum_height minimum_height sampling')
###############################################################################


def getBiomassHeaderOffsetSize(stack_composition):

    return _HEADER_FORMAT_SIZE


def writeBiomassHeader(product_folder, channel_idx, date_str, lon_min, lon_max, lat_min, lat_max, unique_acq_id):

    if not isinstance(date_str, str):
        error_message = 'date_str should be an Utc string, not a PreciseDateTime object'
        logging.error(error_message)
        raise ValueError(error_message)

    if len(date_str) != _date_str_STR_LEN:
        error_message = 'date_str has a different length from an Utc date'
        logging.error(error_message)
        raise ValueError(error_message)

    # encode all the strings with the _STRING_ENCODING format
    encoded_date_str = date_str.encode(_STRING_ENCODING)
    encoded_unique_acq_id = unique_acq_id.encode(_STRING_ENCODING)

    # fill the struct with all data to write
    packed_data = struct.pack(
        _HEADER_FORMAT, encoded_date_str, lon_min, lon_max, lat_min, lat_max, encoded_unique_acq_id
    )

    raster_info_read = (
        product_folder.get_channel(channel_idx).metadata.get_metadata_channels(0).get_element('RasterInfo')
    )

    if raster_info_read.header_offset_bytes != _HEADER_FORMAT_SIZE:
        error_message = 'Incompatible header offset size, please follow this flow: step 1 -> getBiomassHeaderOffsetSize; step 2 -> append_channel to product_folder with offset from step 1; step 3 -> execute this function'
        logging.error(error_message)
        raise ValueError(error_message)

    raster_file = os.path.join(product_folder.pf_dir_path, raster_info_read.file_name)
    with open(raster_file, 'wb') as raster_fid:
        raster_fid.write(packed_data)


def readBiomassHeader(product_folder, channel_idx):

    data_channel_obj = product_folder.get_channel(channel_idx)
    metadata_obj = data_channel_obj.metadata
    metadatachannel_obj = metadata_obj.get_metadata_channels(0)
    ri = metadatachannel_obj.get_element('RasterInfo')
    raster_file = os.path.join(product_folder.pf_dir_path, ri.file_name)

    date_str, lon_min, lon_max, lat_min, lat_max, unique_acq_id = readBiomassHeader_core(raster_file)

    return date_str, lon_min, lon_max, lat_min, lat_max, unique_acq_id


def readBiomassHeader_core(raster_file):

    # get the product folder RasterInfo and retrive the raster_file name

    # read the raster file (just the header)
    with open(raster_file, 'br') as fid:
        packed_data = fid.read(_HEADER_FORMAT_SIZE)

    encoded_date_str, lon_min, lon_max, lat_min, lat_max, encoded_unique_acq_id = struct.unpack(
        _HEADER_FORMAT, packed_data
    )

    date_str = encoded_date_str.decode(_STRING_ENCODING)
    unique_acq_id = encoded_unique_acq_id.decode(_STRING_ENCODING)

    return date_str, lon_min, lon_max, lat_min, lat_max, unique_acq_id


def decode_unique_acquisition_id_string(unique_acquisition_id_string, output_format='numeric'):

    if output_format == 'string':
        global_cycle_idx = unique_acquisition_id_string[3:5]
        heading_deg = unique_acquisition_id_string[8:14]
        rg_swath_idx = unique_acquisition_id_string[20:22]
        rg_sub_swath_idx = unique_acquisition_id_string[30:32]
        az_swath_idx = unique_acquisition_id_string[38:40]
        baseline_idx = unique_acquisition_id_string[45:]

    elif output_format == 'numeric':
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

    return global_cycle_idx, heading_deg, rg_swath_idx, rg_sub_swath_idx, az_swath_idx, baseline_idx


def getBinaryNameFromChannelIDX(pf_name, channel_idx):

    if not type(channel_idx) == type(1) or channel_idx < 1:
        error_message = 'Input channel_idx should be an integer number greater than zero'
        logging.error(error_message)
        raise ValueError(error_message)

    channel_idx_str = str(channel_idx)
    number_of_zeros = 4 - len(channel_idx_str)
    zeros_str = '0' * number_of_zeros
    binary_name = pf_name + '_' + zeros_str + channel_idx_str

    return binary_name


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


def write_chains_input_file(input_params, input_file_xml):
    """Write the input file of the inner chains"""

    chain_id = input_params.chain_id
    InputsL2_Item = Element('InputsL2' + chain_id)
    InputsL2_Item.set('version', _XML_VERSION)

    L1cRepository = SubElement(InputsL2_Item, 'L1cRepository')
    L1cRepository.text = input_params.L1c_repository

    # L1cProductList
    L1cProductList = SubElement(InputsL2_Item, 'L1cProductList')

    for stack_id, stack_composition_list in input_params.stack_composition.items():

        L1cProduct = SubElement(L1cProductList, 'L1cProduct')
        L1cProduct.set('unique_stack_id', stack_id)

        for pf_name in stack_composition_list:

            Acquisition = SubElement(L1cProduct, 'Acquisition')
            Acquisition.text = pf_name

    # AuxiliaryProductList
    AuxiliaryProductList = SubElement(InputsL2_Item, 'AuxiliaryProductList')

    if len(input_params.ECEF_grid_file_names):
        for stack_id, stack_composition_list in input_params.stack_composition.items():

            L1cStack = SubElement(AuxiliaryProductList, 'L1cStack')
            L1cStack.set('unique_stack_id', stack_id)

            Geometry = SubElement(L1cStack, 'Geometry')

            ECEFGrid = SubElement(Geometry, 'ECEFGrid')
            ECEFGrid.text = input_params.ECEF_grid_file_names[stack_id]

            kz = SubElement(Geometry, 'kz')
            kz.text = input_params.kz_file_names[stack_id]

            RadarDistances = SubElement(Geometry, 'RadarDistances')
            RadarDistances.text = input_params.slant_range_distances_file_names[stack_id]

            OffNadirAngle = SubElement(Geometry, 'OffNadirAngles')
            OffNadirAngle.text = input_params.off_nadir_angle_file_names[stack_id]

            Slope = SubElement(Geometry, 'Slope')
            Slope.text = input_params.slope_file_names[stack_id]

            ReferenceHeight = SubElement(Geometry, 'ReferenceHeight')
            ReferenceHeight.text = input_params.reference_height_file_names[stack_id]

            CalibrationScreens = SubElement(L1cStack, 'CalibrationScreens')
            CalibrationScreens.text = input_params.calibration_screens_file_names[stack_id]

    if chain_id == 'FD':
        AverageCovariance = SubElement(AuxiliaryProductList, 'AverageCovarianceFolder')
        AverageCovariance.text = input_params.average_covariance_folder

    DEM = SubElement(AuxiliaryProductList, 'DEM')
    DEM.text = input_params.dem_folder

    if chain_id == 'AGB':
        ReferenceAGB = SubElement(AuxiliaryProductList, 'ReferenceAGB')
        ReferenceAGB.text = input_params.reference_agb_folder

        ForestHeight = SubElement(AuxiliaryProductList, 'ForestHeight')
        ForestHeight.text = input_params.forest_height_folder

    if chain_id == 'FH':
        SystemDecorrelationFunction = SubElement(AuxiliaryProductList, 'SystemDecorrelationFunction')
        SystemDecorrelationFunction.text = input_params.system_decorrelation_fun_folder

    if not 'TOMO' in chain_id:
        ForestMask = SubElement(AuxiliaryProductList, 'ForestMask')
        ForestMask.text = input_params.forest_mask_catalogue_folder

    # OutputFolder
    OutputFolder = SubElement(InputsL2_Item, 'OutputFolder')
    OutputFolder.text = input_params.output_folder

    # GeographicGridSampling
    GeographicGridSampling = SubElement(InputsL2_Item, 'GeographicGridSampling')
    GeographicGridSampling.text = str(input_params.geographic_grid_sampling)
    GeographicGridSampling.set('unit', 'm')

    # write to file
    ElementTree_indent(InputsL2_Item)
    tree = ElementTree(InputsL2_Item)
    tree.write(open(input_file_xml, 'w'), encoding='unicode')


def write_chains_configuration_file(configuration_params, configuration_file_xml):
    """Write the input file of the inner chains"""

    chain_id = configuration_params.chain_id

    ConfigurationL2_Item = Element('ConfigurationL2' + chain_id)
    ConfigurationL2_Item.set('version', _XML_VERSION)

    if chain_id == 'TOMO':
        InterpolateStack = SubElement(ConfigurationL2_Item, 'InterpolateStack')
        RegularBaselineGrid = SubElement(InterpolateStack, 'RegularBaselineGrid')
        RegularBaselineGrid.text = str(configuration_params.interpolate_stack.regular_baseline_grid_flag)

    # Ground Cancellation
    if chain_id == 'AGB' or chain_id == 'FD':

        ground_canc_params = configuration_params.ground_cancellation
        GroundCancellation = SubElement(ConfigurationL2_Item, 'GroundCancellation')
        MultiMaster = SubElement(GroundCancellation, 'MultiMaster')
        MultiMaster.text = str(ground_canc_params.multi_master_flag)
        EnhancedForestHeight = SubElement(GroundCancellation, 'EnhancedForestHeight')
        EnhancedForestHeight.text = str(ground_canc_params.enhanced_forest_height)
        ModelBasedEqualization = SubElement(GroundCancellation, 'ModelBasedEqualization')
        ModelBasedEqualization.text = str(ground_canc_params.equalization_flag)
        if (
            not ModelBasedEqualization.text == '1'
            and not ModelBasedEqualization.text == '2'
            and not ModelBasedEqualization.text == '3'
        ):
            error_str = 'Configuration flag "ModelBasedEqualization" value "{}" not valid, choose among "1", "2" or "3", where "1"->always OFF; "2"->always ON; "3"->ON only if two acquisitions are present'.format(
                ModelBasedEqualization.text
            )
            logging.error(error_str)
            raise

    params_curr = getattr(configuration_params, chain_id)

    if chain_id == 'AGB':

        Estimate_elem = SubElement(ConfigurationL2_Item, 'Estimate' + chain_id)

        NLooks = SubElement(Estimate_elem, 'NLooks')
        NLooks.text = str(params_curr.intermediate_ground_averaging)

        ProductResolution = SubElement(Estimate_elem, 'ProductResolution')
        ProductResolution.text = str(params_curr.product_resolution)

        DistanceROI = SubElement(Estimate_elem, 'DistanceROI')
        DistanceROI.text = str(params_curr.distance_sampling_area)

        BlockROI = SubElement(Estimate_elem, 'BlockROI')
        BlockROI.text = str(params_curr.ROIs_per_block)

        OverlapROI = SubElement(Estimate_elem, 'OverlapROI')
        OverlapROI.text = str(params_curr.ROIs_per_blocks_overlap)

        model_params = params_curr.model_parameters
        ModelParameters = SubElement(Estimate_elem, 'ModelParameters')

        InitValue = SubElement(ModelParameters, 'InitValue')
        A_elem = SubElement(InitValue, 'A')
        A_elem.text = str(model_params.init_value.A)
        alpha_elem = SubElement(InitValue, 'alpha')
        alpha_elem.text = str(model_params.init_value.alpha)
        n_elem = SubElement(InitValue, 'n')
        n_elem.text = str(model_params.init_value.n)

        ChangesAcrossStacks = SubElement(ModelParameters, 'ChangesAcrossStacks')
        A_elem = SubElement(ChangesAcrossStacks, 'A')
        A_elem.text = str(model_params.changes_across_stacks.A)
        alpha_elem = SubElement(ChangesAcrossStacks, 'alpha')
        alpha_elem.text = str(model_params.changes_across_stacks.alpha)
        n_elem = SubElement(ChangesAcrossStacks, 'n')
        n_elem.text = str(model_params.changes_across_stacks.n)

        ChangesAcrossPol = SubElement(ModelParameters, 'ChangesAcrossPol')
        A_elem = SubElement(ChangesAcrossPol, 'A')
        A_elem.text = str(model_params.changes_across_pol.A)
        alpha_elem = SubElement(ChangesAcrossPol, 'alpha')
        alpha_elem.text = str(model_params.changes_across_pol.alpha)
        n_elem = SubElement(ChangesAcrossPol, 'n')
        n_elem.text = str(model_params.changes_across_pol.n)

    if chain_id == 'FH':

        Estimate_elem = SubElement(ConfigurationL2_Item, 'Estimate' + chain_id)

        SpectralShiftFiltering = SubElement(Estimate_elem, 'SpectralShiftFiltering')
        SpectralShiftFiltering.text = str(params_curr.spectral_shift_filtering)

        ProductResolution = SubElement(Estimate_elem, 'ProductResolution')
        ProductResolution.text = str(params_curr.product_resolution)

        model_params = params_curr.model_parameters
        ModelParameters = SubElement(Estimate_elem, 'ModelParameters')

        MaximumHeight = SubElement(ModelParameters, 'MaximumHeight')
        MaximumHeight.text = str(model_params.maximum_height)

        NumberOfExtinctionValue = SubElement(ModelParameters, 'NumberOfExtinctionValue')
        NumberOfExtinctionValue.text = str(model_params.number_of_extinction_value)

        NumberOfGroundVolumeRatioValue = SubElement(ModelParameters, 'NumberOfGroundVolumeRatioValue')
        NumberOfGroundVolumeRatioValue.text = str(model_params.number_of_ground_volume_ratio_value)

        NumberOfTemporalDecorrelationValue = SubElement(ModelParameters, 'NumberOfTemporalDecorrelationValue')
        NumberOfTemporalDecorrelationValue.text = str(model_params.number_of_temporal_decorrelation_value)

        MedianFactor = SubElement(Estimate_elem, 'MedianFactor')
        MedianFactor.text = str(params_curr.median_factor)

    if chain_id == 'FD':

        Estimate_elem = SubElement(ConfigurationL2_Item, 'ChangeDetection')

        ProductResolution = SubElement(Estimate_elem, 'ProductResolution')
        ProductResolution.text = str(params_curr.product_resolution)

        ConfidenceLevel = SubElement(Estimate_elem, 'ConfidenceLevel')
        ConfidenceLevel.text = str(params_curr.confidence_level)

    if chain_id == 'TOMO_FH' or chain_id == 'TOMO':

        ProductResolution = SubElement(ConfigurationL2_Item, 'ProductResolution')
        ProductResolution.text = str(params_curr.product_resolution)

        vertical_range = configuration_params.vertical_range

        VerticalRange = SubElement(ConfigurationL2_Item, 'VerticalRange')

        MaximumHeight = SubElement(VerticalRange, 'MaximumHeight')
        MaximumHeight.text = str(vertical_range.maximum_height)

        MinimumHeight = SubElement(VerticalRange, 'MinimumHeight')
        MinimumHeight.text = str(vertical_range.minimum_height)

        sampling = SubElement(VerticalRange, 'sampling')
        sampling.text = str(vertical_range.sampling)

    if chain_id == 'TOMO_FH':

        Estimate_elem = SubElement(ConfigurationL2_Item, 'EstimateFH')

        EnableSuperResolution = SubElement(Estimate_elem, 'EnableSuperResolution')
        EnableSuperResolution.text = str(params_curr.enable_super_resolution)

        RegularizationNoiseFactor = SubElement(Estimate_elem, 'RegularizationNoiseFactor')
        RegularizationNoiseFactor.text = str(params_curr.regularization_noise_factor)

        PowerThreshold = SubElement(Estimate_elem, 'PowerThreshold')
        PowerThreshold.text = str(params_curr.power_threshold)

        MedianFactor = SubElement(Estimate_elem, 'MedianFactor')
        MedianFactor.text = str(params_curr.median_factor)

    SaveBreakpoints = SubElement(ConfigurationL2_Item, 'SaveBreakpoints')
    SaveBreakpoints.text = str(configuration_params.save_breakpoints)

    # write to file
    ElementTree_indent(ConfigurationL2_Item)
    tree = ElementTree(ConfigurationL2_Item)
    tree.write(open(configuration_file_xml, 'w'), encoding='unicode')


def parse_biomassL2_main_input_file(input_file_xml):
    """Parse the input file, configuration file and current campaign params file
    and store all the needed information in the dataset data_stack class"""

    # 1) Parse input_file_xml -------------------------------------------------
    tree = ET.parse(input_file_xml)
    root = tree.getroot()

    L1c_repository = root.find('L1cRepository').text
    if not os.path.exists(L1c_repository):
        error_message = 'Main input file: the specified L1c Repository folder does not exist '
        logging.error(error_message)
        raise ValueError(error_message)

    L1c_aux_data_repository = root.find('AuxiliaryProductsFolder').text
    if not os.path.exists(L1c_aux_data_repository):
        error_message = 'Main input file: the specified AuxiliaryProductsFolder does not exist: '
        logging.error(error_message)
        raise ValueError(error_message)

    dem_folder_found = False
    aux_list = os.listdir(L1c_aux_data_repository)
    for aux_name in aux_list:
        if aux_name == 'DEM':
            dem_folder_found = True
    if not dem_folder_found:
        error_message = (
            'Main input file: the specified AuxiliaryProductsFolder should contain AT LEAST the sub-folder "DEM"'
        )
        logging.error(error_message)
        raise ValueError(error_message)

    output_folder = root.find('OutputFolder').text

    try:
        geographic_grid_sampling = float(root.find('GeographicGridSampling').text)
    except:
        error_message = 'Main input file: GeographicGridSampling should be a numeric value'
        logging.error(error_message)
        raise ValueError(error_message)
    if geographic_grid_sampling <= 0:
        error_message = 'Main input file: GeographicGridSampling should be a positive numeric value [m]'
        logging.error(error_message)
        raise ValueError(error_message)

    L2Product_Item = root.find('L2Product')
    proc_flags_AGB = bool_from_string(L2Product_Item.find('AGB').text)
    proc_flags_FD = bool_from_string(L2Product_Item.find('FD').text)
    proc_flags_FH = bool_from_string(L2Product_Item.find('FH').text)
    proc_flags_TOMO = bool_from_string(L2Product_Item.find('TOMO').text)
    proc_flags_TOMO_FH = bool_from_string(L2Product_Item.find('TOMO_FH').text)

    proc_flags_struct = proc_flags(proc_flags_AGB, proc_flags_FH, proc_flags_TOMO_FH, proc_flags_FD, proc_flags_TOMO)

    L1cDates = root.findall('L1cDate')
    if len(L1cDates) != 2:
        error_message = 'Main input file: there should be #2 L1cDates (one with value="start" and one with value="stop"'
        logging.error(error_message)
        raise ValueError(error_message)
    for date_UTC_Item in L1cDates:

        L1cDate = date_UTC_Item.text
        attrib = date_UTC_Item.attrib['value']

        if 'start' in attrib:
            L1cDate_start = PreciseDateTime().set_from_utc_string(L1cDate)
        elif 'stop' in attrib:
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
    GeographicBoundary_Item = root.find('GeographicBoundary')
    Points = GeographicBoundary_Item.findall('Point')
    if len(Points) < 3:
        error_message = 'Main input file: at least #3 Geograplic Boundary points should be specified'
        logging.error(error_message)
        raise ValueError(error_message)

    for Point_Item in Points:

        Latitude_str = Point_Item.find('Latitude').text
        Longitude_str = Point_Item.find('Longitude').text
        Lat = float(Latitude_str)
        Lon = float(Longitude_str)
        if Lat < -90 or Lat > 90:
            error_message = 'Main input file: Latitude should be a number from -90 to +90 [deg] '
            logging.error(error_message)
            raise ValueError(error_message)
        if Lon < -180 or Lon > 180:
            error_message = 'Main input file: Longitude should be a number from -180 to +180 [deg] '
            logging.error(error_message)
            raise ValueError(error_message)

        point_dict = {'Latitude': Lat, 'Longitude': Lon}

        point_list.append(point_dict)
        del point_dict

    geographic_boundary_struct = geographic_boundary(point_list)

    stacks_to_find_list = []
    Stacks_to_find_Item = root.find('Stacks_to_find')
    if Stacks_to_find_Item != None:
        stacks_to_find = Stacks_to_find_Item.findall('Stack_id')
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


def parse_fd_lookup_table(lookup_table_file_name_xml):
    """Parse the look up tables uded to indicize the fnf mask and covariance at each FD step"""

    lut_dict = {}
    # 1) Parse input_file_xml -------------------------------------------------
    tree = ET.parse(lookup_table_file_name_xml)
    root = tree.getroot()

    L1cStacks = root.findall('L1cStack')
    for L1cStack in L1cStacks:

        lookup_table_struct = namedtuple('lookup_table', 'boundaries dates file_names')
        stack_id = L1cStack.attrib['stack_id']

        # extract all the bondaries coordinates from current stack (one date for stack):
        geographic_boundaries_curr = geographic_boundaries(None, None, None, None)

        GeographicBoundaries_Item = L1cStack.find('GeographicBoundaries')
        geographic_boundaries_curr.lat_min = GeographicBoundaries_Item.find('latMin').text
        geographic_boundaries_curr.lat_max = GeographicBoundaries_Item.find('latMax').text
        geographic_boundaries_curr.lon_min = GeographicBoundaries_Item.find('lonMin').text
        geographic_boundaries_curr.lon_max = GeographicBoundaries_Item.find('lonMax').text

        # extract all the dates from current stack (one date for each step):
        dates_dict = {}
        fnames_dict = {}
        try:
            if L1cStack.getchildren()[0].tag == 'GeographicBoundaries':
                table_type = L1cStack.getchildren()[1].tag
            else:
                table_type = L1cStack.getchildren()[0].tag
        except:
            raise ImportError(' table xml is not valid')

        if not table_type == 'fnf' and not table_type == 'covariance':
            raise ImportError(' table xml is not valid, element "' + table_type + '" not recognized')

        fd_steps = L1cStack.findall(table_type)
        for step_curr in fd_steps:
            step_number = step_curr.attrib['step']
            dates_dict[step_number] = step_curr.find('Date').text
            fnames_dict[step_number] = step_curr.find('FileName').text

        # fill the optupt dictionary with the current stack boundaries and dates:
        lookup_table_struct.boundaries = geographic_boundaries_curr
        lookup_table_struct.dates = dates_dict
        lookup_table_struct.file_names = fnames_dict
        lut_dict[stack_id] = lookup_table_struct

        del lookup_table_struct, dates_dict

    return lut_dict


def write_fd_lookup_table(lookup_table_file_name_xml, table_type, lut_dict):
    """Write the lookup table for fd chain fnf and covariance  outputs"""

    # table_type: string 'fnf' or 'covariance'
    if not table_type == 'fnf' and not table_type == 'covariance':
        raise ValueError(
            '#3th iinput should be a string of value "fnf" or "covariance": "' + table_type + '" is not recognized'
        )

    if table_type == 'fnf':
        lut_Item = Element('FDlookUpTableFNFmasks')
    elif table_type == 'covariance':
        lut_Item = Element('FDlookUpTableCovariances')

    lut_Item.set('version', _XML_VERSION)

    for stack_id in lut_dict.keys():

        lookup_table_struct = lut_dict[stack_id]

        L1cStack = SubElement(lut_Item, 'L1cStack')
        L1cStack.set('stack_id', stack_id)

        GeographicBoundaries = SubElement(L1cStack, 'GeographicBoundaries')
        GeographicBoundaries.set('unit', 'deg')

        lat_min = SubElement(GeographicBoundaries, 'latMin')
        lat_min.text = lookup_table_struct.boundaries.lat_min
        lat_max = SubElement(GeographicBoundaries, 'latMax')
        lat_max.text = lookup_table_struct.boundaries.lat_max
        lon_min = SubElement(GeographicBoundaries, 'lonMin')
        lon_min.text = lookup_table_struct.boundaries.lon_min
        lon_max = SubElement(GeographicBoundaries, 'lonMax')
        lon_max.text = lookup_table_struct.boundaries.lon_max

        for step_str_idx in lookup_table_struct.dates.keys():

            dates_and_fnames_element = SubElement(L1cStack, table_type)
            dates_and_fnames_element.set('step', step_str_idx)

            date = SubElement(dates_and_fnames_element, 'Date')
            date.text = lookup_table_struct.dates[step_str_idx]
            date.set('unit', 'Utc')

            file_name = SubElement(dates_and_fnames_element, 'FileName')
            file_name.text = lookup_table_struct.file_names[step_str_idx]

    # write to file
    ElementTree_indent(lut_Item)
    tree = ElementTree(lut_Item)
    tree.write(open(lookup_table_file_name_xml, 'w'), encoding='unicode')


def bool_from_string(in_string):

    in_string = in_string[0].upper() + in_string[1:].lower()

    if 'True' in in_string or 'False' in in_string:
        bool_value = ast.literal_eval(in_string)

    elif '0' in in_string or '1' in in_string:
        bool_value = bool(ast.literal_eval(in_string))
    else:
        error_message = 'Main input file: the L2product AGB, FH, FD and TOMO values should be boolean: True, False, 0 or 1 are supported'
        logging.error(error_message)
        raise ValueError(error_message)

    return bool_value


def parse_chains_input_file(input_file_xml):
    """Parse the input file of inner chains(AGB, FH, FD, TOMO FH and TOMO)"""

    # 1) Parse input_file_xml -------------------------------------------------
    tree = ET.parse(input_file_xml)
    root = tree.getroot()
    chain_id = root.tag[8:]
    if chain_id != 'AGB' and chain_id != 'FH' and chain_id != 'FD' and chain_id != 'TOMO_FH' and chain_id != 'TOMO':
        error_message = 'The provided input file is not an inner biomassL2 chain input (AGB, FH, FD, TOMO FH or TOMO)'
        logging.error(error_message)
        raise ValueError(error_message)

    L1cRepository_Item = root.find('L1cRepository')
    L1c_repository = L1cRepository_Item.text

    L1cProductList_Item = root.find('L1cProductList')
    L1cProduct_Items = L1cProductList_Item.findall('L1cProduct')
    number_of_stacks = len(L1cProduct_Items)
    stack_composition = {}
    for L1cProduct_Item in L1cProduct_Items:

        stack_id = L1cProduct_Item.attrib['unique_stack_id']

        Acquisition_Items = L1cProduct_Item.findall('Acquisition')
        Acquisitions_list = []
        for Acquisition_Item in Acquisition_Items:
            Acquisitions_list.append(Acquisition_Item.text)

        stack_composition[stack_id] = Acquisitions_list

    AuxiliaryProductList_Item = root.find('AuxiliaryProductList')
    DEM_Item = AuxiliaryProductList_Item.find('DEM')
    dem_folder = DEM_Item.text

    if chain_id == 'FD':
        AVG_Item = AuxiliaryProductList_Item.find('AverageCovarianceFolder')
        average_covariance_folder = AVG_Item.text
    else:
        average_covariance_folder = ''

    L1cStack_Items = AuxiliaryProductList_Item.findall('L1cStack')
    if len(L1cStack_Items) != number_of_stacks:
        error_message = (
            chain_id
            + ' Input File: the number of L1cStack elements should be equal to the number of L1cProduct elements'
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

        stack_id = L1cStack_Item.attrib['unique_stack_id']
        if not stack_id in stack_composition.keys():
            error_message = (
                chain_id
                + ' Input File: the "stack_id" values of L1cStack elelemts does not match the "stack_id" values of L1cProduct elements'
            )
            logging.error(error_message)
            raise ValueError(error_message)

        Geometry_Item = L1cStack_Item.find('Geometry')

        ECEF_grid_file_names[stack_id] = Geometry_Item.find('ECEFGrid').text
        kz_file_names[stack_id] = Geometry_Item.find('kz').text
        slant_range_distances_file_names[stack_id] = Geometry_Item.find('RadarDistances').text
        off_nadir_angle_file_names[stack_id] = Geometry_Item.find('OffNadirAngles').text
        slope_file_names[stack_id] = Geometry_Item.find('Slope').text
        reference_height_file_names[stack_id] = Geometry_Item.find('ReferenceHeight').text

        if L1cStack_Item.find('CalibrationScreens') != None:
            calibration_screens_file_names[stack_id] = L1cStack_Item.find('CalibrationScreens').text
        elif chain_id == 'TOMO' or chain_id == 'TOMO_FH':
            error_message = chain_id + ' Input File: CalibrationScreens are mandatory for ' + chain_id + ' chain.'
            logging.error(error_message)
            raise ValueError(error_message)

    reference_agb_folder = ''
    system_decorrelation_fun_folder = ''
    forest_mask_catalogue_folder = ''
    forest_height_folder = ''

    if chain_id == 'AGB' and AuxiliaryProductList_Item.find('ReferenceAGB').text:
        reference_agb_folder = AuxiliaryProductList_Item.find('ReferenceAGB').text
    if chain_id == 'AGB' and AuxiliaryProductList_Item.find('ForestHeight').text:
        forest_height_folder = AuxiliaryProductList_Item.find('ForestHeight').text
    if chain_id == 'FH':
        if AuxiliaryProductList_Item.find('SystemDecorrelationFunction'):
            system_decorrelation_fun_folder = AuxiliaryProductList_Item.find('SystemDecorrelationFunction').text
        else:
            system_decorrelation_fun_folder = ''

    if not 'TOMO' in chain_id:
        forest_mask_catalogue_folder = AuxiliaryProductList_Item.find('ForestMask').text

    output_folder = root.find('OutputFolder').text

    geographic_grid_sampling = float(root.find('GeographicGridSampling').text)

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


def check_chains_input_file(proc_inputs):

    if not os.path.exists(proc_inputs.L1c_repository):
        error_message = ' Input file L1cRepository does not exist: ' + proc_inputs.L1c_repository
        logging.error(error_message)
        raise RuntimeError(error_message)

    for key in proc_inputs.stack_composition.keys():
        for pf_name in proc_inputs.stack_composition[key]:
            fullPath = os.path.join(proc_inputs.L1c_repository, pf_name)
            if not os.path.exists(fullPath):
                error_message = ' Input file Acquisition ' + pf_name + ' does not exist: ' + fullPath
                logging.error(error_message)
                raise RuntimeError(error_message)

    if not os.path.exists(proc_inputs.dem_folder):
        error_message = ' Input file AuxiliaryProductList DEM does not exist: ' + proc_inputs.dem_folder
        logging.error(error_message)
        raise RuntimeError(error_message)

    for ecef_name in proc_inputs.ECEF_grid_file_names:
        if not os.path.exists(ecef_name):
            error_message = ' Input file AuxiliaryProductList ECEFgrid does not exist: ' + ecef_name
            logging.error(error_message)
            raise RuntimeError(error_message)

    for kz_name in proc_inputs.kz_file_names:
        if not os.path.exists(kz_name):
            error_message = ' Input file AuxiliaryProductList KZ does not exist: ' + kz_name
            logging.error(error_message)
            raise RuntimeError(error_message)

    for inc_name in proc_inputs.off_nadir_angle_file_names:
        if not os.path.exists(inc_name):
            error_message = ' Input file AuxiliaryProductList OffNadirAngle does not exist: ' + inc_name
            logging.error(error_message)
            raise RuntimeError(error_message)

    for h_name in proc_inputs.reference_height_file_names:
        if not os.path.exists(h_name):
            error_message = 'Input file AuxiliaryProductList ReferenceHeight does not exist: ' + h_name
            logging.error(error_message)
            raise RuntimeError(error_message)

    chain_id = proc_inputs.chain_id
    if chain_id == 'AGB':
        if not os.path.exists(proc_inputs.reference_agb_folder):
            error_message = (
                ' Input file AuxiliaryProductList ReferenceAGB does not exist: ' + proc_inputs.reference_agb_folder
            )
            logging.error(error_message)
            raise RuntimeError(error_message)
    if chain_id == 'FD':
        if not os.path.exists(proc_inputs.average_covariance_folder):
            error_message = (
                ' Input  AuxiliaryProductList AverageCovarianceFolder does not exist: '
                + proc_inputs.average_covariance_folder
            )
            logging.error(error_message)
            raise RuntimeError(error_message)
    if chain_id == 'FH':
        if not os.path.exists(proc_inputs.system_decorrelation_fun_folder):
            error_message = (
                ' Input file AuxiliaryProductList SystemDecorrelationFunction does not exist: '
                + proc_inputs.system_decorrelation_fun_folder
            )
            logging.error(error_message)
            raise RuntimeError(error_message)

    if not os.path.exists(proc_inputs.forest_mask_catalogue_folder):
        error_message = (
            ' Input file AuxiliaryProductList forest non forest mask catalogue folder does not exist: '
            + proc_inputs.forest_mask_catalogue_folder
        )
        logging.error(error_message)
        raise RuntimeError(error_message)

    if len(proc_inputs.geographic_grid_sampling) == 0 or proc_inputs.geographic_grid_sampling <= 0:
        error_message = (
            ' Input file GeographicGridSampling should be a positive number: ' + proc_inputs.geographic_grid_sampling
        )
        logging.error(error_message)
        raise RuntimeError(error_message)

    if not os.path.exists(proc_inputs.forest_mask_catalogue_folder):
        error_message = (
            ' Input file AuxiliaryProductList ForestMask does not exist: ' + proc_inputs.forest_mask_catalogue_folder
        )
        logging.error(error_message)
        raise RuntimeError(error_message)


def parse_chains_configuration_file(configuration_file_xml):
    """Parse the chains configuration files for AGB, FD, FH, TOMO, TOMO_FH"""

    tree = ET.parse(configuration_file_xml)
    root = tree.getroot()
    chain_id = root.tag[15:]

    # interpolate_stack
    if chain_id == 'TOMO':
        InterpolateStack_Item = root.find('InterpolateStack')
        regular_baseline_grid_flag = bool_from_string(InterpolateStack_Item.find('RegularBaselineGrid').text)
        interpolate_stack = interpolate_params(regular_baseline_grid_flag)
    else:
        interpolate_stack = None

    # ground_cancellation
    if chain_id == 'AGB' or chain_id == 'FD':
        GroundCancellation_Item = root.find('GroundCancellation')
        multi_master_flag = bool_from_string(GroundCancellation_Item.find('MultiMaster').text)
        enhanced_forest_height = float(GroundCancellation_Item.find('EnhancedForestHeight').text)
        equalization_flag = GroundCancellation_Item.find('ModelBasedEqualization').text
        if not equalization_flag == '1' and not equalization_flag == '2' and not equalization_flag == '3':
            error_str = 'Configuration flag "ModelBasedEqualization" value "{}" not valid, choose among "1", "2" or "3", where "1"->always OFF; "2"->always ON; "3"->ON only if two acquisitions are present'.format(
                equalization_flag
            )
            logging.error(error_str)
            raise
        ground_cancellation = ground_canc_params(multi_master_flag, enhanced_forest_height, equalization_flag)
    else:
        ground_cancellation = None

    # vertical_range
    if chain_id == 'TOMO' or chain_id == 'TOMO_FH':
        VerticalRange_Item = root.find('VerticalRange')
        maximum_height = float(VerticalRange_Item.find('MaximumHeight').text)
        minimum_height = float(VerticalRange_Item.find('MinimumHeight').text)
        sampling = float(VerticalRange_Item.find('sampling').text)
        vertical_range = vertical_range_params(maximum_height, minimum_height, sampling)
    else:
        vertical_range = None

    enable_resampling = bool_from_string(root.find('EnableResampling').text)
    compute_geometry = bool_from_string(root.find('ComputeGeometry').text)
    apply_calibration_screen = bool_from_string(root.find('ApplyCalibrationScreen').text)
    DEM_flattening = bool_from_string(root.find('DEMflattening').text)
    multilook_heading_correction = bool_from_string(root.find('MultilookHeadingCorrection').text)
    save_breakpoints = bool_from_string(root.find('SaveBreakpoints').text)
    delete_temporary_files = bool_from_string(root.find('DeleteTemporaryFiles').text)

    if chain_id == 'AGB':
        chain_field_name = 'EstimateAGB'
    elif chain_id == 'FD':
        chain_field_name = 'ChangeDetection'
    elif chain_id == 'FH' or chain_id == 'TOMO_FH':
        chain_field_name = 'EstimateFH'
    elif chain_id == 'TOMO':
        chain_field_name = ''

    chain_field_Item = root.find(chain_field_name)

    if chain_id == 'TOMO' or chain_id == 'TOMO_FH':
        product_resolution = float(root.find('ProductResolution').text)
    elif chain_id == 'AGB':
        product_resolution = float(root.find(chain_field_name).find('product_resolution').text)
    else:
        product_resolution = float(chain_field_Item.find('ProductResolution').text)

    AGB = None
    FH = None
    FD = None
    TOMO_FH = None
    TOMO = None
    if chain_id == 'AGB':

        number_of_tests = float(chain_field_Item.find('number_of_tests').text)
        intermediate_ground_averaging = float(chain_field_Item.find('intermediate_ground_averaging').text)

        distance_sampling_area = float(chain_field_Item.find('distance_sampling_area').text)
        fraction_of_roi_per_test = int(chain_field_Item.find('fraction_of_roi_per_test').text)
        fraction_of_cal_per_test = int(chain_field_Item.find('fraction_of_cal_per_test').text)

        parameter_block_size = float(chain_field_Item.find('parameter_block_size').text)
        distance_parameter_block = float(chain_field_Item.find('distance_parameter_block').text)
        min_number_of_rois = int(chain_field_Item.find('min_number_of_rois').text)
        min_number_of_rois_per_stack = int(chain_field_Item.find('min_number_of_rois_per_stack').text)
        min_number_of_cals_per_test = int(chain_field_Item.find('min_number_of_cals_per_test').text)
        min_number_of_rois_per_test = int(chain_field_Item.find('min_number_of_rois_per_test').text)

        model_parameter_Item = chain_field_Item.find('ModelParameters')

        parameter_ranges_Item = model_parameter_Item.find('parameter_ranges')
        agb_model_scaling_limits = [
            float(parameter_ranges_Item.find('agb_model_scaling').find('min').text),
            float(parameter_ranges_Item.find('agb_model_scaling').find('max').text),
        ]
        agb_model_exponent_limits = [
            float(parameter_ranges_Item.find('agb_model_exponent').find('min').text),
            float(parameter_ranges_Item.find('agb_model_exponent').find('max').text),
        ]
        agb_cosine_exponent_limits = [
            float(parameter_ranges_Item.find('agb_cosine_exponent').find('min').text),
            float(parameter_ranges_Item.find('agb_cosine_exponent').find('max').text),
        ]
        agb_estimation_limits = [
            float(parameter_ranges_Item.find('agb_estimation').find('min').text),
            float(parameter_ranges_Item.find('agb_estimation').find('max').text),
        ]

        changes_across_pol_Item = model_parameter_Item.find('ChangesAcrossPolarisation')
        agb_model_scaling = bool_from_string(changes_across_pol_Item.find('agb_model_scaling').text)
        agb_model_exponent = bool_from_string(changes_across_pol_Item.find('agb_model_exponent').text)
        agb_cosine_exponent = bool_from_string(changes_across_pol_Item.find('agb_cosine_exponent').text)
        changes_across_pol = triplet_params(agb_model_scaling, agb_model_exponent, agb_cosine_exponent)

        changes_across_stacks_Item = model_parameter_Item.find('ChangesAcrossStack')
        agb_model_scaling = bool_from_string(changes_across_stacks_Item.find('agb_model_scaling').text)
        agb_model_exponent = bool_from_string(changes_across_stacks_Item.find('agb_model_exponent').text)
        agb_cosine_exponent = bool_from_string(changes_across_stacks_Item.find('agb_cosine_exponent').text)
        changes_across_stacks = triplet_params(agb_model_scaling, agb_model_exponent, agb_cosine_exponent)

        changes_across_global_cycle_Item = model_parameter_Item.find('ChangesAcrossGlobalCycle')
        agb_model_scaling = bool_from_string(changes_across_global_cycle_Item.find('agb_model_scaling').text)
        agb_model_exponent = bool_from_string(changes_across_global_cycle_Item.find('agb_model_exponent').text)
        agb_cosine_exponent = bool_from_string(changes_across_global_cycle_Item.find('agb_cosine_exponent').text)
        changes_across_global_cycle = triplet_params(agb_model_scaling, agb_model_exponent, agb_cosine_exponent)

        changes_across_heading_Item = model_parameter_Item.find('ChangesAcrossHeading')
        agb_model_scaling = bool_from_string(changes_across_heading_Item.find('agb_model_scaling').text)
        agb_model_exponent = bool_from_string(changes_across_heading_Item.find('agb_model_exponent').text)
        agb_cosine_exponent = bool_from_string(changes_across_heading_Item.find('agb_cosine_exponent').text)
        changes_across_heading = triplet_params(agb_model_scaling, agb_model_exponent, agb_cosine_exponent)

        changes_across_swath_Item = model_parameter_Item.find('ChangesAcrossSwath')
        agb_model_scaling = bool_from_string(changes_across_swath_Item.find('agb_model_scaling').text)
        agb_model_exponent = bool_from_string(changes_across_swath_Item.find('agb_model_exponent').text)
        agb_cosine_exponent = bool_from_string(changes_across_swath_Item.find('agb_cosine_exponent').text)
        changes_across_swath = triplet_params(agb_model_scaling, agb_model_exponent, agb_cosine_exponent)

        changes_across_sub_swath_Item = model_parameter_Item.find('ChangesAcrossSubswath')
        agb_model_scaling = bool_from_string(changes_across_sub_swath_Item.find('agb_model_scaling').text)
        agb_model_exponent = bool_from_string(changes_across_sub_swath_Item.find('agb_model_exponent').text)
        agb_cosine_exponent = bool_from_string(changes_across_sub_swath_Item.find('agb_cosine_exponent').text)
        changes_across_subswath = triplet_params(agb_model_scaling, agb_model_exponent, agb_cosine_exponent)

        changes_across_azimuth_Item = model_parameter_Item.find('ChangesAcrossAzimuth')
        agb_model_scaling = bool_from_string(changes_across_azimuth_Item.find('agb_model_scaling').text)
        agb_model_exponent = bool_from_string(changes_across_azimuth_Item.find('agb_model_exponent').text)
        agb_cosine_exponent = bool_from_string(changes_across_azimuth_Item.find('agb_cosine_exponent').text)
        changes_across_azimuth = triplet_params(agb_model_scaling, agb_model_exponent, agb_cosine_exponent)

        agb_model_scaling = bool_from_string(changes_across_stacks_Item.find('agb_model_scaling').text)
        agb_model_exponent = bool_from_string(changes_across_stacks_Item.find('agb_model_exponent').text)
        agb_cosine_exponent = bool_from_string(changes_across_stacks_Item.find('agb_cosine_exponent').text)
        changes_across_stacks = triplet_params(agb_model_scaling, agb_model_exponent, agb_cosine_exponent)

        estimation_valid_values_limits = [
            float(chain_field_Item.find('EstimationValidValuesLimits').find('min').text),
            float(chain_field_Item.find('EstimationValidValuesLimits').find('max').text),
        ]

        model_parameters = agb_model_params(
            agb_model_scaling_limits,
            agb_model_exponent_limits,
            agb_cosine_exponent_limits,
            agb_estimation_limits,
            changes_across_pol,
            changes_across_stacks,
            changes_across_global_cycle,
            changes_across_heading,
            changes_across_swath,
            changes_across_subswath,
            changes_across_azimuth,
        )

        AGB = agb_est_params(
            number_of_tests,
            fraction_of_roi_per_test,
            fraction_of_cal_per_test,
            intermediate_ground_averaging,
            product_resolution,
            distance_sampling_area,
            parameter_block_size,
            distance_parameter_block,
            min_number_of_rois,
            min_number_of_rois_per_stack,
            min_number_of_cals_per_test,
            min_number_of_rois_per_test,
            model_parameters,
            estimation_valid_values_limits,
        )

    elif chain_id == 'FD':

        confidence_level = float(chain_field_Item.find('ConfidenceLevel').text)
        FD = FD_est_params(product_resolution, confidence_level)

    elif chain_id == 'FH':

        estimation_valid_values_limits = [
            float(chain_field_Item.find('EstimationValidValuesLimits').find('min').text),
            float(chain_field_Item.find('EstimationValidValuesLimits').find('max').text),
        ]

        spectral_shift_filtering = bool_from_string(chain_field_Item.find('SpectralShiftFiltering').text)

        kz_thresholds = [
            float(chain_field_Item.find('KZ_thresholds').find('min').text),
            float(chain_field_Item.find('KZ_thresholds').find('max').text),
        ]

        median_factor_temp = float(chain_field_Item.find('MedianFactor').text)
        median_factor = int(median_factor_temp)
        if median_factor_temp - median_factor > 0:
            logging.warning(
                'FH configuration file MedianFactor should be integer, the value that will be used is: {}'.format(
                    median_factor
                )
            )

        model_parameters_Item = chain_field_Item.find('ModelParameters')
        maximum_height = int(model_parameters_Item.find('MaximumHeight').text)
        number_of_extinction_value = int(model_parameters_Item.find('NumberOfExtinctionValue').text)
        number_of_ground_volume_ratio_value = int(model_parameters_Item.find('NumberOfGroundVolumeRatioValue').text)
        number_of_temporal_decorrelation_value = int(
            model_parameters_Item.find('NumberOfTemporalDecorrelationValue').text
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

        FH = FH_est_params(spectral_shift_filtering, product_resolution, kz_thresholds, model_parameters, median_factor)

    elif chain_id == 'TOMO_FH':

        estimation_valid_values_limits = [
            float(chain_field_Item.find('EstimationValidValuesLimits').find('min').text),
            float(chain_field_Item.find('EstimationValidValuesLimits').find('max').text),
        ]

        enable_super_resolution = bool_from_string(chain_field_Item.find('EnableSuperResolution').text)
        regularization_noise_factor = float(chain_field_Item.find('RegularizationNoiseFactor').text)
        power_threshold = float(chain_field_Item.find('PowerThreshold').text)
        median_factor = int(chain_field_Item.find('MedianFactor').text)
        TOMO_FH = TOMO_FH_est_params(
            estimation_valid_values_limits,
            product_resolution,
            enable_super_resolution,
            regularization_noise_factor,
            power_threshold,
            median_factor,
        )

    elif chain_id == 'TOMO':

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


def parse_campaign_params_simple(campaign_params_file_xml, tag, campaign_tag):
    """Parse the campaign params file xml"""
    # Parse each campaign _param.xml file ---------------------------------------
    tree = ET.parse(campaign_params_file_xml)
    root = tree.getroot()

    scene_node = root.find('Scenes')
    carrierFrequency = float(scene_node.attrib['frequency']) * 1000000

    temp = [sc.attrib['pixel_spacing'] for sc in scene_node.findall('Scene') if sc.attrib['name'] in tag]
    pixel_spacing_slant_rg_m = float(temp[0])
    del temp

    if campaign_tag in ['tropisar']:
        pixel_spacing_az_m = 1.245
    elif campaign_tag in ['afrisar_onera']:
        pixel_spacing_az_m = 1.2
    elif campaign_tag in ['afrisar_dlr']:
        pixel_spacing_az_m = 0.90919
    else:
        pixel_spacing_az_m = 1

    temp = [sc.attrib['SLR_start'] for sc in scene_node.findall('Scene') if sc.attrib['name'] in tag]
    SLR_start_m = float(temp[0])
    del temp

    date = [sc.attrib['date'] for sc in scene_node.findall('Scene') if sc.attrib['name'] in tag]

    months = ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DIC']
    month_str = months[int(date[0][3:5]) - 1]

    date_UTC = (
        date[0][:2]
        + '-'
        + month_str
        + '-'
        + date[0][6:10]
        + ' '
        + date[0][11:13]
        + ':'
        + date[0][14:17]
        + ':00.000000000000'
    )
    '23-MAY-1998 08:25:12.123456789012345'

    temp = [sc.attrib['master'] for sc in scene_node.findall('Scene') if sc.attrib['name'] in tag]
    master = temp[0]

    temp = [sc.attrib['heading'] for sc in scene_node.findall('Scene') if sc.attrib['name'] in tag]
    heading = float(temp[0])

    return pixel_spacing_az_m, pixel_spacing_slant_rg_m, date_UTC, SLR_start_m, master, carrierFrequency, heading
