import os
import logging
import pyproj
import shutil
import collections
import numpy as np
from equi7grid.equi7grid import Equi7Grid
from equi7grid.image2equi7grid import image2equi7grid
from skimage.filters.rank import majority as majority_filter
from scipy.interpolate import griddata, interp1d
from scipy.constants import c as LIGHTSPEED
from arepytools.timing.precisedatetime import PreciseDateTime
from shapely.geometry import MultiPoint
from osgeo import gdal, osr
from gdalconst import GA_ReadOnly
from datetime import datetime

from biomassL2.IO_interfaces import (
    readBiomassHeader_core,
    getBinaryNameFromChannelIDX,
    input_params,
    write_chains_input_file,
    raster_info,
    geographic_boundaries,
    decode_unique_acquisition_id_string,
)
from biomassL2.data_import import (
    read_data,
    read_ecef_grid,
    read_auxiliary_multi_channels,
    read_auxiliary_single_channel,
    data_oversample,
    get_data_time_stamp,
)
from biomassL2.utility_functions import epsg_in_to_epsg_out, tandemx_fnf_read
from biomassL2.constants import (
    MAX_EARTH_RADIUS_METERS,
    OVERSAMPLING_FACTOR,
    EPSG_CODE_ECEF,
    EPSG_CODE_LLA,
    SATELLITE_VELOCITY,
)
from biomassL2 import geometric_lib

from arepytools.io.productfolder import ProductFolder
from arepytools.io.metadata import DataSetInfo
from arepytools.io.metadata import RasterInfo
from arepytools.geometry.generalsarorbit import create_general_sar_orbit


def start_logging(output_folder, proc_flags, log_level):
    # CRITICAL 50
    # ERROR 40
    # WARNING 30
    # INFO 20
    # DEBUG 10
    # NOTSET 0

    if log_level == 'DEBUG':
        level_to_set = logging.DEBUG
    elif log_level == 'INFO':
        level_to_set = logging.INFO
    elif log_level == 'WARNING':
        level_to_set = logging.WARNING
    elif log_level == 'ERROR':
        level_to_set = logging.ERROR

    log_file_name = os.path.join(output_folder, 'biomassL2.log')

    logging.basicConfig(
        handlers=[logging.FileHandler(log_file_name, mode='w', encoding='utf-8'), logging.StreamHandler()],
        level=level_to_set,
        format='%(asctime)s - %(levelname)s | %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
    )

    logging.getLogger('matplotlib.font_manager').disabled = True

    logging.info(' --BIOMASS L2 Processor-- ')
    logging.info('Following chains will be executed: ')
    if proc_flags.AGB:
        logging.info('AGB')
    if proc_flags.FD:
        logging.info('FD')
    if proc_flags.FH:
        logging.info('FH (Interferometric phase)')
    if proc_flags.TOMO_FH:
        logging.info('FH (Tomographic phase)')
    if proc_flags.TOMO:
        logging.info('TOMO (Tomographic cube generation)')
    logging.info(' \n')

    return log_file_name


def get_auxiliary_paths(stack_composition, L1c_aux_data_repository):
    # get the names of the folders containing of the auxiliary data: the names
    # shouls match the default ones

    ECEFGRID_folder = os.path.join(L1c_aux_data_repository, 'Geometry', 'ECEFGRID')
    KZ_folder = os.path.join(L1c_aux_data_repository, 'Geometry', 'KZ')
    off_nadir_angle_folder = os.path.join(L1c_aux_data_repository, 'Geometry', 'OffNadirAngles')
    reference_height_folder = os.path.join(L1c_aux_data_repository, 'Geometry', 'ReferenceHeight')
    slant_range_distances = os.path.join(L1c_aux_data_repository, 'Geometry', 'SlantRangeDistances')
    slope_folder = os.path.join(L1c_aux_data_repository, 'Geometry', 'Slope')

    forest_mask_catalogue_folder = os.path.join(L1c_aux_data_repository, 'ForestMask')
    dem_folder = os.path.join(L1c_aux_data_repository, 'DEM')
    reference_agb_folder = os.path.join(L1c_aux_data_repository, 'ReferenceAGB')
    average_covariance_folder = os.path.join(L1c_aux_data_repository, 'AverageCovariance')
    system_decorrelation_fun_folder = os.path.join(L1c_aux_data_repository, 'SystemDecorrelationFunction')
    calibration_screens_folder = os.path.join(L1c_aux_data_repository, 'CalibrationScreens')
    forest_height_folder = os.path.join(L1c_aux_data_repository, 'ForestHeight')

    # Those are stack based auxiliary:
    ECEF_grid_file_names = {}
    kz_file_names = {}
    off_nadir_angle_file_names = {}
    reference_height_file_names = {}
    slant_range_distances_file_names = {}
    slope_file_names = {}
    calibration_screens_file_names = {}

    for stack_id in stack_composition.keys():

        ECEF_grid_file_names[stack_id] = os.path.join(r'' + ECEFGRID_folder, stack_id)
        kz_file_names[stack_id] = os.path.join(r'' + KZ_folder, stack_id)
        off_nadir_angle_file_names[stack_id] = os.path.join(r'' + off_nadir_angle_folder, stack_id)
        reference_height_file_names[stack_id] = os.path.join(r'' + reference_height_folder, stack_id)
        slant_range_distances_file_names[stack_id] = os.path.join(r'' + slant_range_distances, stack_id)
        slope_file_names[stack_id] = os.path.join(r'' + slope_folder, stack_id)
        calibration_screens_file_names[stack_id] = os.path.join(r'' + calibration_screens_folder, stack_id)

    return (
        ECEF_grid_file_names,
        kz_file_names,
        off_nadir_angle_file_names,
        reference_height_file_names,
        slant_range_distances_file_names,
        slope_file_names,
        forest_mask_catalogue_folder,
        dem_folder,
        reference_agb_folder,
        average_covariance_folder,
        system_decorrelation_fun_folder,
        calibration_screens_folder,
        forest_height_folder,
        calibration_screens_file_names,
    )


def data_select_by_date_and_boundaries(main_input_struct):

    geographic_boundaries_curr = geographic_boundaries(None, None, None, None)

    # from all the database of data, selectonly the ones which matches the
    # input parameters and return a dictionary with the stacks to be used

    start_date = main_input_struct.L1c_date.start
    stop_date = main_input_struct.L1c_date.stop

    # the stack_composition will contain Stack_IDs and Stack Compositions for the available dataset at desired dates and boundaries
    stacks_composition = {}
    # stack_composition_dates_match and stack_composition_boundaries_match are used only for warning display
    # in case no data is found in the desired dates+boudaries, but just for one of the two
    stacks_composition_only_dates_match = {}
    stacks_composition_only_boundaries_match = {}
    stacks_composition_not_enough_acq = {}
    stacks_geographic_boundaries = {}

    # initialize the list which will contain all the header_biomass extracted realizations

    # cycle over all the data set and read the binary header: we read just first channel, no need for all of them
    channel_idx = 1
    pf_names_list = os.listdir(main_input_struct.L1c_repository)

    for pf_idx, pf_name in enumerate(pf_names_list):
        date_is_matching = False
        boundaries_are_matching = False
        stack_is_matching = False
        num_acq_are_enough = False

        if len(pf_name) != 47:
            error_msg = 'Unrecognized uniqie acquisition ID (it will be skipped): "{}" from repository {}; the correct format is "GC_XX_H_XXX.XX_RGSW_XX_RGSBSW_X_AZSW_XX_BSL_XX" where each "X" is an integer'.format(
                pf_name, main_input_struct.L1c_repository
            )
            logging.error(error_msg)
            continue

        output_pf_full_name = os.path.join(main_input_struct.L1c_repository, pf_name)

        bin_file_name = getBinaryNameFromChannelIDX(pf_name, channel_idx)
        raster_file = os.path.join(output_pf_full_name, bin_file_name)
        if not os.path.exists(raster_file):
            logging.warning(
                'The folder "'
                + bin_file_name
                + '" inside L1C repository is NOT a product folder and it will be skipped \n'
            )
            continue

        # read the current header

        datestr, lon_min, lon_max, lat_min, lat_max, unique_acq_id = readBiomassHeader_core(raster_file)

        uniqie_stack_id = unique_acq_id[0:40]

        # first test:check the date:
        date = PreciseDateTime().set_from_utc_string(datestr)

        if pf_idx == 0:
            edge_values = {
                'date_min': date,
                'date_max': date,
                'lon_min': lon_min,
                'lon_max': lon_max,
                'lat_min': lat_min,
                'lat_max': lat_max,
            }

        else:
            edge_values['date_min'] = min(edge_values['date_min'], date)
            edge_values['date_max'] = max(edge_values['date_max'], date)
            edge_values['lon_min'] = min(edge_values['lon_min'], lon_min)
            edge_values['lon_max'] = max(edge_values['lon_max'], lon_max)
            edge_values['lat_min'] = min(edge_values['lat_min'], lat_min)
            edge_values['lat_max'] = max(edge_values['lat_max'], lat_max)

        if uniqie_stack_id in stacks_geographic_boundaries.keys():

            stacks_geographic_boundaries[uniqie_stack_id].lon_min = min(
                stacks_geographic_boundaries[uniqie_stack_id].lon_min, lon_min
            )
            stacks_geographic_boundaries[uniqie_stack_id].lon_max = max(
                stacks_geographic_boundaries[uniqie_stack_id].lon_max, lon_max
            )
            stacks_geographic_boundaries[uniqie_stack_id].lat_min = min(
                stacks_geographic_boundaries[uniqie_stack_id].lat_min, lat_min
            )
            stacks_geographic_boundaries[uniqie_stack_id].lat_max = max(
                stacks_geographic_boundaries[uniqie_stack_id].lat_max, lat_max
            )
        else:
            geographic_boundaries_per_stack = geographic_boundaries(None, None, None, None)
            geographic_boundaries_per_stack.lon_min = lon_min
            geographic_boundaries_per_stack.lon_max = lon_max
            geographic_boundaries_per_stack.lat_min = lat_min
            geographic_boundaries_per_stack.lat_max = lat_max
            stacks_geographic_boundaries[uniqie_stack_id] = geographic_boundaries_per_stack

        # check 2/4 search for data inside the user dates boundaries:
        date_is_matching = check_on_dates(start_date, stop_date, date)

        # check 3/4 search for data inside the user ground coordinates boundaries:
        boundaries_are_matching = check_on_boundaries(
            lon_min, lon_max, lat_min, lat_max, main_input_struct.geographic_boundary.point
        )

        # check 4/4 (optional) ckeck if data matches one of the user specified stack:
        stack_is_matching = check_on_stacks_to_find(uniqie_stack_id, main_input_struct.stacks_to_find)

        if date_is_matching and boundaries_are_matching and stack_is_matching:  # and num_acq_are_enough:
            # store information for data that satisfy all previous checks:
            if not uniqie_stack_id in stacks_composition.keys():
                stacks_composition[uniqie_stack_id] = [unique_acq_id]
            else:
                stacks_composition[uniqie_stack_id].append(unique_acq_id)

            try:
                geographic_boundaries_curr.lon_min = min(geographic_boundaries_curr.lon_min, lon_min)
                geographic_boundaries_curr.lon_max = max(geographic_boundaries_curr.lon_max, lon_max)
                geographic_boundaries_curr.lat_min = min(geographic_boundaries_curr.lat_min, lat_min)
                geographic_boundaries_curr.lat_max = max(geographic_boundaries_curr.lat_max, lat_max)
            except:
                geographic_boundaries_curr.lon_min = lon_min
                geographic_boundaries_curr.lon_max = lon_max
                geographic_boundaries_curr.lat_min = lat_min
                geographic_boundaries_curr.lat_max = lat_max

        elif date_is_matching and stack_is_matching and num_acq_are_enough:
            # store information for warning display
            if not uniqie_stack_id in stacks_composition_only_dates_match.keys():
                stacks_composition_only_dates_match[uniqie_stack_id] = [unique_acq_id]
            else:
                stacks_composition_only_dates_match[uniqie_stack_id].append(unique_acq_id)

        elif boundaries_are_matching and stack_is_matching and num_acq_are_enough:
            # store information for warning display
            if not uniqie_stack_id in stacks_composition_only_boundaries_match.keys():
                stacks_composition_only_boundaries_match[uniqie_stack_id] = [unique_acq_id]
            else:
                stacks_composition_only_boundaries_match[uniqie_stack_id].append(unique_acq_id)

    # return the headers_list and some useful messages to the user
    num_stacks = len(stacks_composition)
    if num_stacks == 0:
        logging.warning('\n')
        logging.warning('Cannot find any valid data:')
        if len(stacks_composition_only_dates_match) or len(stacks_composition_only_boundaries_match):

            logging.warning('    be aware that data In L1cRepository are confined inside following boundaries:')
            logging.warning('    minimun date is       {}'.format(edge_values['date_min']))
            logging.warning('    maximun date is       {}'.format(edge_values['date_max']))
            logging.warning('    minimun Longitude is  {}'.format(edge_values['lon_min']))
            logging.warning('    maximum Longitude is  {}'.format(edge_values['lon_max']))
            logging.warning('    minimun Latitude  is  {}'.format(edge_values['lat_min']))
            logging.warning('    maximum Latitude  is  {} \n'.format(edge_values['lat_max']))
            logging.warning('    Instead the user asked for following boundaries:')
            logging.warning('    user min date         {}'.format(start_date))
            logging.warning('    user max date         {}'.format(stop_date))

            point_dict = main_input_struct.geographic_boundary.point

            logging.warning('    user min Longitude    {}'.format(min([x['Longitude'] for x in point_dict])))
            logging.warning('    user max Longitude    {}'.format(max([x['Longitude'] for x in point_dict])))
            logging.warning('    user min Latitude     {}'.format(min([x['Latitude'] for x in point_dict])))
            logging.warning('    user max Latitude     {} \n'.format(max([x['Latitude'] for x in point_dict])))

            if len(stacks_composition_only_dates_match):
                logging.warning(
                    '    Following stacks are confined in the desired L1cDates but not in the desired geographic boundaries'
                )
                logging.warning('    {} \n'.format(list(stacks_composition_only_dates_match.keys())))

            if len(stacks_composition_only_boundaries_match):
                logging.warning(
                    '    Following stacks are confined in the desired geographic boundaries but not in the desired L1cDates'
                )
                logging.warning('    {} \n'.format(list(stacks_composition_only_boundaries_match.keys())))

        if len(stacks_composition_not_enough_acq):
            logging.warning(
                'Following stacks contain just ONE acquisition (a minimum of TWO are needed to have ONE baseline)'
            )
            logging.warning('{} \n'.format(list(stacks_composition_not_enough_acq.keys())))

        if (
            not len(stacks_composition_only_dates_match)
            and not len(stacks_composition_only_boundaries_match)
            and not len(stacks_composition_not_enough_acq)
        ):
            logging.warning('Cannot find any data, the L1cRepository is empty')
            pf_idx = 0

    else:
        logging.info('Found #{} stacks: {}'.format(num_stacks, list(stacks_composition.keys())))
        logging.info('Stacks composition: \n')
        for stack_curr, pf_list in stacks_composition.items():
            pf_list.sort()
            logging.info('Stack ' + stack_curr + ' composition ( unique acquisitions id ):')
            for pf_name in pf_list:
                logging.info('    ' + pf_name)

    stacks_composition = collections.OrderedDict(sorted(stacks_composition.items()))

    return stacks_composition, geographic_boundaries_curr, stacks_geographic_boundaries


def check_on_dates(start_date, stop_date, date):
    # core function to check the dates match

    if date >= start_date and date <= stop_date:
        date_is_matching = True
    else:
        date_is_matching = False

    return date_is_matching


def check_on_boundaries(lon_min, lon_max, lat_min, lat_max, point_dict):
    # core function to check if the data falls inside the user geographic boundaries

    # create the geographic polygon containing the user input searching area:
    data_longitudes = [x['Longitude'] for x in point_dict]
    data_latitudes = [x['Latitude'] for x in point_dict]
    data_lat_lon_coords = tuple(zip(data_latitudes, data_longitudes))
    data_polygon_obj = MultiPoint(data_lat_lon_coords).convex_hull

    # create the geographic polygon containing the current data
    user_lat_lon_coords = tuple(zip([lat_min, lat_min, lat_max, lat_max], [lon_min, lon_max, lon_min, lon_max]))
    user_polygon_obj = MultiPoint(user_lat_lon_coords).convex_hull

    # check if the data polygon is at least partially intersected by the user polygon
    boundaries_are_matching = user_polygon_obj.intersects(data_polygon_obj)

    return boundaries_are_matching


def check_on_stacks_to_find(uniqie_stack_id, stacks_to_find):
    # core function: this is optional, in the sense that stacks_to_find can be empty, in this case don't check on stacks (set it to True)
    if len(stacks_to_find):
        if uniqie_stack_id in stacks_to_find:
            stack_is_matching = True
        else:
            stack_is_matching = False
    else:
        stack_is_matching = True

    return stack_is_matching


def write_chains_input_file_main(output_folder, main_input_struct, stack_composition):
    # write the input files of each inner chain

    (
        ECEF_grid_file_names,
        kz_file_names,
        off_nadir_angle_file_names,
        reference_height_file_names,
        slant_range_distances_file_names,
        slope_file_names,
        forest_mask_catalogue_folder,
        dem_folder,
        reference_agb_folder,
        average_covariance_folder,
        system_decorrelation_fun_folder,
        calibration_screens_folder,
        forest_height_folder,
        calibration_screens_file_names,
    ) = get_auxiliary_paths(stack_composition, main_input_struct.L1c_aux_data_repository)

    if main_input_struct.proc_flags.AGB:

        output_folder_sub_chain = os.path.join(output_folder, 'AGB')
        os.makedirs(output_folder_sub_chain)
        AGB_input_file_xml = os.path.join(output_folder_sub_chain, 'InputFile.xml')

        input_params_obj = input_params(
            'AGB',
            main_input_struct.L1c_repository,
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
            output_folder_sub_chain,
            main_input_struct.geographic_grid_sampling,
            None,  # system_decorrelation_fun_folder
            None,  # average_covariance_folder
            forest_height_folder,
        )

        write_chains_input_file(input_params_obj, AGB_input_file_xml)
    else:
        AGB_input_file_xml = ''

    if main_input_struct.proc_flags.FD:

        output_folder_sub_chain = os.path.join(output_folder, 'FD')
        os.makedirs(output_folder_sub_chain)
        FD_input_file_xml = os.path.join(output_folder_sub_chain, 'InputFile.xml')

        input_params_obj = input_params(
            'FD',
            main_input_struct.L1c_repository,
            stack_composition,
            ECEF_grid_file_names,
            kz_file_names,
            slant_range_distances_file_names,
            off_nadir_angle_file_names,
            slope_file_names,
            reference_height_file_names,
            calibration_screens_file_names,
            dem_folder,
            None,  # reference_agb_folder
            forest_mask_catalogue_folder,
            output_folder_sub_chain,
            main_input_struct.geographic_grid_sampling,
            None,  # system_decorrelation_fun_folder
            average_covariance_folder,
            None,
        )  # forest_height_folder

        write_chains_input_file(input_params_obj, FD_input_file_xml)
    else:
        FD_input_file_xml = ''

    if main_input_struct.proc_flags.FH:

        output_folder_sub_chain = os.path.join(output_folder, 'FH')
        os.makedirs(output_folder_sub_chain)
        FH_input_file_xml = os.path.join(output_folder_sub_chain, 'InputFile.xml')

        input_params_obj = input_params(
            'FH',
            main_input_struct.L1c_repository,
            stack_composition,
            ECEF_grid_file_names,
            kz_file_names,
            slant_range_distances_file_names,
            off_nadir_angle_file_names,
            slope_file_names,
            reference_height_file_names,
            calibration_screens_file_names,
            dem_folder,
            None,  # reference_agb_folder
            forest_mask_catalogue_folder,
            output_folder_sub_chain,
            main_input_struct.geographic_grid_sampling,
            system_decorrelation_fun_folder,
            None,  # average_covariance_folder
            None,
        )  # forest_height_folder

        write_chains_input_file(input_params_obj, FH_input_file_xml)
    else:
        FH_input_file_xml = ''

    if main_input_struct.proc_flags.TOMO_FH:

        output_folder_sub_chain = os.path.join(output_folder, 'TOMO_FH')
        os.makedirs(output_folder_sub_chain)

        stack_composition_TOMO = {}
        for stack_id, acquisitions_list in stack_composition.items():
            if len(acquisitions_list) >= 3:
                stack_composition_TOMO[stack_id] = stack_composition[stack_id]

        len_stacks_all = stack_composition.keys()
        len_stacks_tomo = stack_composition_TOMO.keys()
        if len_stacks_tomo:
            if len_stacks_tomo < len_stacks_all:
                logging.info(
                    'The stacks with less than #3 acquisitions cannot be used to perform TOMO FH estimation \
                and they will be removed: #{} stacks among #{} totaL, have been removed'.format(
                        len_stacks_tomo, len_stacks_all
                    )
                )

            TOMO_FH_input_file_xml = os.path.join(output_folder_sub_chain, 'InputFile.xml')

            input_params_obj = input_params(
                'TOMO_FH',
                main_input_struct.L1c_repository,
                stack_composition,
                ECEF_grid_file_names,
                kz_file_names,
                slant_range_distances_file_names,
                off_nadir_angle_file_names,
                slope_file_names,
                reference_height_file_names,
                calibration_screens_file_names,
                dem_folder,
                None,  # reference_agb_folder
                None,  # forest_mask_catalogue_folder
                output_folder_sub_chain,
                main_input_struct.geographic_grid_sampling,
                None,  # system_decorrelation_fun_folder
                None,  # average_covariance_folder
                None,
            )  # forest_height_folder

            write_chains_input_file(input_params_obj, TOMO_FH_input_file_xml)
        else:
            TOMO_FH_input_file_xml = ''
            logging.warning('cannot find any stack with more than #2 acquisitions: TOMO FH estimation will be disabled')
            main_input_struct.proc_flags.TOMO_FH = False
    else:
        TOMO_FH_input_file_xml = ''

    if main_input_struct.proc_flags.TOMO:

        output_folder_sub_chain = os.path.join(output_folder, 'TOMO')
        os.makedirs(output_folder_sub_chain)

        stack_composition_TOMO = {}
        for stack_id, acquisitions_list in stack_composition.items():
            if len(acquisitions_list) >= 3:
                stack_composition_TOMO[stack_id] = stack_composition[stack_id]

        len_stacks_all = stack_composition.keys()
        len_stacks_tomo = stack_composition_TOMO.keys()
        if len_stacks_tomo:
            if len_stacks_tomo < len_stacks_all:
                logging.info(
                    'The stacks with less than #3 acquisitions cannot be used to perform TOMO FH estimation \
                and they will be removed: #{} stacks among #{} totaL, have been removed'.format(
                        len_stacks_tomo, len_stacks_all
                    )
                )

            TOMO_input_file_xml = os.path.join(output_folder_sub_chain, 'InputFile.xml')

            input_params_obj = input_params(
                'TOMO',
                main_input_struct.L1c_repository,
                stack_composition,
                ECEF_grid_file_names,
                kz_file_names,
                slant_range_distances_file_names,
                off_nadir_angle_file_names,
                slope_file_names,
                reference_height_file_names,
                calibration_screens_file_names,
                dem_folder,
                None,  # reference_agb_folder
                None,  # forest_mask_catalogue_folder
                output_folder_sub_chain,
                main_input_struct.geographic_grid_sampling,
                None,  # system_decorrelation_fun_folder
                None,  # average_covariance_folder
                None,
            )  # forest_height_folder

            write_chains_input_file(input_params_obj, TOMO_input_file_xml)

        else:
            TOMO_input_file_xml = ''
            logging.warning('cannot find any stack with more than #2 acquisitions: TOMO estimation will be disabled')
            main_input_struct.proc_flags.TOMO = False

    else:
        TOMO_input_file_xml = ''

    return AGB_input_file_xml, FD_input_file_xml, FH_input_file_xml, TOMO_FH_input_file_xml, TOMO_input_file_xml


def choose_equi7_sampling(product_resolution, geographic_grid_sampling):

    if geographic_grid_sampling <= product_resolution / 2:
        equi7_sampling = geographic_grid_sampling
    else:
        equi7_sampling = product_resolution / 2
        warning_message = 'user rquested Geographic Grid Sampling is {} [m]; it cannot be greater than (product_resolution)/2 = {} [m] \n'.format(
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


def check_if_path_exists(path, file_folder_str):
    if not os.path.exists(path):
        log_error_str = file_folder_str + ' ' + path + ' does not exist'

        logging.error(log_error_str)
        raise Exception(log_error_str)


def read_and_oversample_data(L1c_repository, acquisitions_pf_names, enable_resampling):
    # this function calls the read_data followed by the data_oversample
    # (which oversamples only when needed and if enabled )

    beta0_calibrated = {}
    for pf_name in acquisitions_pf_names:
        logging.info('    loading ' + pf_name + '...')
        # every data in the current stack has same spacing, same resolution and same number of samples and lines
        (
            beta0_calibrated[pf_name],
            num_samples,
            num_lines,
            pixel_spacing_slant_rg,
            pixel_spacing_az,
            resolution_m_slant_rg,
            resolution_m_az,
            master_id,
            lines_start_utc,
        ) = read_data(L1c_repository, pf_name)

    ### ALL chains: oversampling data:
    # needed whenever computing covariance or detected data (as notch for AGB)
    # check

    # backup needed for the resampling in auxiliary data
    raster_info_orig = raster_info(
        num_samples,
        num_lines,
        pixel_spacing_slant_rg,
        pixel_spacing_az,
        resolution_m_slant_rg,
        resolution_m_az,
        lines_start_utc,
    )

    if enable_resampling:

        beta0_calibrated, num_samples, pixel_spacing_slant_rg, num_lines, pixel_spacing_az = data_oversample(
            beta0_calibrated, OVERSAMPLING_FACTOR, raster_info_orig
        )

        logging.info('all data loaded.\n')

    # filling output structure:
    raster_info_os = raster_info(
        num_samples,
        num_lines,
        pixel_spacing_slant_rg,
        pixel_spacing_az,
        resolution_m_slant_rg,
        resolution_m_az,
        lines_start_utc,
    )

    return beta0_calibrated, master_id, raster_info_os, raster_info_orig


def check_if_geometry_auxiliaries_are_present(
    file_names,
    stack_id,
    acquisitions_pf_names,
    read_ecef=True,
    read_off_nadir=True,
    read_slope=True,
    read_kz=True,
    read_ref_h=True,
    read_dist=True,
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
            aux_are_present_flag = False

    if read_slope:
        pf_name = os.path.basename(file_names.slope_file_names[stack_id])
        folder_name = os.path.dirname(file_names.slope_file_names[stack_id])
        if not os.path.exists(os.path.join(folder_name, pf_name)):
            aux_are_present_flag = False

    return aux_are_present_flag


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
        logging.info('Loading auxiliary data: ECEFGRID...')
        pf_name = os.path.basename(file_names.ECEF_grid_file_names[stack_id])
        folder_name = os.path.dirname(file_names.ECEF_grid_file_names[stack_id])
        ecef_grid = read_ecef_grid(folder_name, pf_name)

        if not ecef_grid is None and enable_resampling:
            ecef_grid = data_oversample(ecef_grid, OVERSAMPLING_FACTOR, raster_info)[0]
        if not ecef_grid is None:
            logging.info('...ECEFGRID read.\n')
        else:
            logging.warning(
                'Since "ECEFGRID/'
                + stack_id
                + '" folder is missing into AuxiliaryProductsFolder, geometric library needs to be called.\n'
            )

    if read_kz:
        logging.info('Loading auxiliary data: KZ...')
        pf_name = os.path.basename(file_names.kz_file_names[stack_id])
        folder_name = os.path.dirname(file_names.kz_file_names[stack_id])
        kz = read_auxiliary_multi_channels(folder_name, pf_name, acquisitions_pf_names)

        if not kz is None and enable_resampling:
            kz = data_oversample(kz, OVERSAMPLING_FACTOR, raster_info)[0]
        if not kz is None:
            logging.info('...KZ read.\n')
        else:
            logging.warning(
                'Since "KZ/'
                + stack_id
                + '" folder is missing into AuxiliaryProductsFolder, geometric library needs to be called.\n'
            )

    if read_off_nadir:
        logging.info('Loading auxiliary data: off nadir angles...')
        pf_name = os.path.basename(file_names.off_nadir_angle_file_names[stack_id])
        folder_name = os.path.dirname(file_names.off_nadir_angle_file_names[stack_id])
        off_nadir_angle_rad = read_auxiliary_multi_channels(folder_name, pf_name)

        if not off_nadir_angle_rad is None and enable_resampling:
            off_nadir_angle_rad = data_oversample(off_nadir_angle_rad, OVERSAMPLING_FACTOR, raster_info)[0]
        if not off_nadir_angle_rad is None:
            logging.info('...off nadir angles read.\n')
        else:
            logging.warning(
                'Since "OffNadirAngles/'
                + stack_id
                + '" folder is missing into AuxiliaryProductsFolder, geometric library needs to be called.\n'
            )

    if read_ref_h:
        logging.info(
            'Loading auxiliary data: reference height...'
        )  # not needed, use just for comparison with the estimation
        pf_name = os.path.basename(file_names.reference_height_file_names[stack_id])
        folder_name = os.path.dirname(file_names.reference_height_file_names[stack_id])
        reference_height = read_auxiliary_single_channel(
            folder_name, pf_name
        )  # it is the dtm in slant_range-azimuth reference

        if not reference_height is None and enable_resampling:
            reference_height = data_oversample(reference_height, OVERSAMPLING_FACTOR, raster_info)[0]
        if not reference_height is None:
            logging.info('...reference height read.\n')
        else:
            logging.warning(
                'Since "ReferenceHeight/'
                + stack_id
                + '" folder is missing into AuxiliaryProductsFolder, geometric library needs to be called.\n'
            )

    if read_dist:
        logging.info('Loading auxiliary data: slant range distances...')
        pf_name = os.path.basename(file_names.slant_range_distances_file_names[stack_id])
        folder_name = os.path.dirname(file_names.slant_range_distances_file_names[stack_id])
        R = read_auxiliary_multi_channels(folder_name, pf_name)

        if not R is None and enable_resampling:
            R = data_oversample(R, OVERSAMPLING_FACTOR, raster_info)[0]
        if not R is None:
            logging.info('...slant range distances read.\n')
        else:
            logging.warning(
                'Since "SlantRangeDistances/'
                + stack_id
                + '" folder is missing into AuxiliaryProductsFolder, geometric library needs to be called.\n'
            )

    if read_slope:
        logging.info('Loading auxiliary data: slope...')
        pf_name = os.path.basename(file_names.slope_file_names[stack_id])
        folder_name = os.path.dirname(file_names.slope_file_names[stack_id])
        slope = read_auxiliary_single_channel(folder_name, pf_name)

        if not slope.any() is None and enable_resampling:
            slope = data_oversample(slope, OVERSAMPLING_FACTOR, raster_info)[0]
        if not slope is None:
            logging.info('...slope read.\n')
        else:
            logging.warning(
                'Since "Slopes/'
                + stack_id
                + '" folder is missing into AuxiliaryProductsFolder, geometric library needs to be called.\n'
            )

    if read_average_cov:
        logging.info('Loading auxiliary data: Average covariance matrix...')
        pf_name = os.path.basename(file_names.average_covariance_folder)
        folder = os.path.dirname(file_names.average_covariance_folder)
        average_covariance = read_auxiliary_multi_channels(folder, pf_name, acquisitions_pf_names)

        if not average_covariance is None and enable_resampling:
            average_covariance = data_oversample(average_covariance, OVERSAMPLING_FACTOR, raster_info)[0]
        if not average_covariance is None:
            logging.info('... Average covariance matrix read.\n')

    if read_cal_screens:
        logging.info('Loading auxiliary data: Calibration Screens...')
        pf_name = os.path.basename(file_names.calibration_screens_file_names[stack_id])
        folder = os.path.dirname(file_names.calibration_screens_file_names[stack_id])
        calibration_screens, cal_screens_raster_info = read_auxiliary_multi_channels(
            folder, pf_name, acquisitions_pf_names, read_raster_info=True
        )

        if not calibration_screens is None and enable_resampling:
            calibration_screens = data_oversample(calibration_screens, OVERSAMPLING_FACTOR, raster_info)[0]
        if not calibration_screens is None:
            logging.info('... Calibration Screens read.\n')

    if read_FH:
        logging.info('Loading auxiliary data: forest estimated height..')
        pf_name = os.path.basename(file_names.forest_height_folder)
        folder = os.path.dirname(file_names.forest_height_folder)
        forest_height = read_auxiliary_single_channel(folder, pf_name, acquisitions_pf_names)

        if not forest_height is None and enable_resampling:
            forest_height = data_oversample(forest_height, OVERSAMPLING_FACTOR, raster_info)[0]
        logging.info('...done.\n')
        if not forest_height is None:
            logging.info('...forest estimated height read.\n')

    if read_reference_agb:
        logging.info('Loading auxiliary data: reference agb..')
        pf_name = os.path.basename(file_names.reference_agb_folder)
        folder = os.path.dirname(file_names.reference_agb_folder)
        reference_agb = read_auxiliary_single_channel(folder, pf_name, acquisitions_pf_names)

        if not reference_agb is None and enable_resampling:
            reference_agb = data_oversample(reference_agb, OVERSAMPLING_FACTOR, raster_info)[0]
        if not reference_agb is None:
            logging.info('...reference agb read.\n')

    if read_sys_dec_fun:
        logging.info('Loading auxiliary data: system decorrelation function...')
        pf_name = os.path.basename(file_names.system_decorrelation_fun_folder)
        folder = os.path.dirname(file_names.system_decorrelation_fun_folder)
        system_decorr_fun = read_auxiliary_single_channel(folder, pf_name, acquisitions_pf_names)

        if not system_decorr_fun is None and enable_resampling:
            system_decorr_fun = data_oversample(system_decorr_fun, OVERSAMPLING_FACTOR, raster_info)[0]
        if not system_decorr_fun is None:
            logging.info('...system decorrelation function read.\n')

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


def evaluate_estimation_quality_matrix(data_in_shape):
    # Placemark for the quality estimation to be defined

    logging.warning(
        'The quality estimation of the product is still to be defined: zeros will be placed in the quality layer '
    )
    quality_matrix = np.zeros(data_in_shape)

    return quality_matrix


def geocoding_init(ecef_grid, rg_vec_subs, az_vec_subs, min_spacing_m):

    ### Convert ECEFGRID pixel to pixel from X+Y+Z to L+L+H (with for example epsg_in_to_epsg_out): the result will NOT be a uniform grid of lat lon
    logging.info('convert ECEFGRID pixel to pixel from X Y Z to Lat Lon Alt...')
    lon_in, lat_in, alt_in = epsg_in_to_epsg_out(
        ecef_grid['X'][rg_vec_subs, :][:, az_vec_subs],
        ecef_grid['Y'][rg_vec_subs, :][:, az_vec_subs],
        ecef_grid['Z'][rg_vec_subs, :][:, az_vec_subs],
        EPSG_CODE_ECEF,
        EPSG_CODE_LLA,
    )

    lat_lon_step = np.rad2deg(min_spacing_m / MAX_EARTH_RADIUS_METERS)

    logging.info('    geocoding latitude and longitude step used: {} [deg]'.format(lat_lon_step))
    lat_lon_gap = 0  # 1*lat_lon_step

    # latitude should be decrescent
    used_lat_step = -np.abs(lat_lon_step)
    lat_regular_vector = np.arange(np.nanmax(lat_in) - lat_lon_gap, np.nanmin(lat_in) + lat_lon_gap, used_lat_step)
    # longitude should be crescent
    used_lon_step = np.abs(lat_lon_step)
    lon_regular_vector = np.arange(np.nanmin(lon_in) - lat_lon_gap, np.nanmax(lon_in) + lat_lon_gap, used_lon_step)

    lonMeshed_out, latMeshed_out = np.meshgrid(lon_regular_vector, lat_regular_vector)

    valid_values_mask = np.invert(np.isnan(lat_in) + np.isnan(lon_in) + np.isnan(alt_in))

    return (
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
    )


def geocoding(data_to_geocode, lon_in, lat_in, lonMeshed_out, latMeshed_out, valid_values_mask, interp_method='linear'):
    f = griddata(
        np.array([lat_in[valid_values_mask].flatten(), lon_in[valid_values_mask].flatten()]).T,
        data_to_geocode[valid_values_mask].flatten(),
        np.array([latMeshed_out.flatten(), lonMeshed_out.flatten()]).T,
        method=interp_method,
    )  # interpolate to fill holes
    data_ground = f.reshape(latMeshed_out.shape)

    return data_ground


def tiff_formatter(
    data_in, out_fname, geotransform, gdal_data_format, projection=None, multi_layers_tiff=False, time_tag=None
):

    if '.tiff' in out_fname:
        out_fname = out_fname[0:-5]
    elif '.tif' in out_fname:
        out_fname = out_fname[0:-4]

    if isinstance(data_in, list) and multi_layers_tiff:
        # write multi layer data in same tiff

        if isinstance(geotransform[0], list):
            geotransform = geotransform[0]

        out_tiff_fname = out_fname + '.tif'
        num_layers = len(data_in)
        # formats and saves the input data in GEO-TIFF
        if type(data_in[0]) == str:
            data_temp = np.load(data_in[0])
            Nx, Ny = data_temp.shape
            del data_temp
        else:
            Nx, Ny = data_in[0].shape

        driver = gdal.GetDriverByName("GTiff")

        outdata = driver.Create(out_tiff_fname, Ny, Nx, num_layers, gdal_data_format)
        if time_tag:
            outdata.SetMetadata({'time_tag': time_tag}, "TIFFTAG_DATETIME")

        if projection:
            outdata.SetProjection(projection)

        else:

            srs = osr.SpatialReference()
            srs.ImportFromEPSG(np.int(EPSG_CODE_LLA[5:]))
            outdata.SetProjection(srs.ExportToWkt())

        outdata.SetGeoTransform(geotransform)

        for idx, data in enumerate(data_in):

            if type(data) == str:
                outdata.GetRasterBand(idx + 1).WriteArray(np.load(data))

            else:
                outdata.GetRasterBand(idx + 1).WriteArray(data)

        outdata.FlushCache()  ##saves to disk!!
        outdata = None

    elif isinstance(data_in, list) and not multi_layers_tiff:
        # write each data in a different tiff

        out_tiff_fname = []
        for idx, data in enumerate(data_in):
            out_tiff_fname.append(out_fname + '_fnftile' + str(idx) + '.tif')

            # formats and saves the input data in GEO-TIFF
            Nx, Ny = data.shape

            driver = gdal.GetDriverByName("GTiff")

            outdata = driver.Create(out_tiff_fname[idx], Ny, Nx, 1, gdal_data_format)
            if time_tag:
                outdata.SetMetadata({'time_tag': time_tag}, "TIFFTAG_DATETIME")

            if projection:
                outdata.SetProjection(projection)

            else:
                srs = osr.SpatialReference()
                srs.ImportFromEPSG(np.int(EPSG_CODE_LLA[5:]))
                outdata.SetProjection(srs.ExportToWkt())

            outdata.SetGeoTransform(geotransform[idx])
            outdata.GetRasterBand(1).WriteArray(data)
            outdata.FlushCache()  ##saves to disk!!
            outdata = None

    else:

        # write the single input data to tiff
        out_tiff_fname = out_fname + '.tif'

        # formats and saves the input data in GEO-TIFF
        Nx, Ny = data_in.shape

        driver = gdal.GetDriverByName("GTiff")

        outdata = driver.Create(out_tiff_fname, Ny, Nx, 1, gdal_data_format)
        if time_tag:
            outdata.SetMetadata({'time_tag': time_tag}, "TIFFTAG_DATETIME")

        if projection:
            outdata.SetProjection(projection)

        else:
            srs = osr.SpatialReference()
            srs.ImportFromEPSG(np.int(EPSG_CODE_LLA[5:]))
            outdata.SetProjection(srs.ExportToWkt())

        outdata.SetGeoTransform(geotransform)
        outdata.GetRasterBand(1).WriteArray(data_in)
        outdata.FlushCache()  ##saves to disk!!
        outdata = None

    return out_tiff_fname


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

        # remember: geotransform = [ lon_start, lon_step, 0, lat_start, 0, lat_step]
        lon0 = geotransform_list[idx_fnf][0] + geotransform_list[idx_fnf][1] * len_x / 2
        lon1 = lon0 + geotransform_list[idx_fnf][1]
        lat0 = geotransform_list[idx_fnf][3] + geotransform_list[idx_fnf][5] * len_x / 2
        lat1 = lat0 + geotransform_list[idx_fnf][5]

        geod = pyproj.Geod(ellps='WGS84')

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
        error_msg = (
            'Input FNF Mask (TANDEM-X) date of "{}" is bigger than input stack data minimum date of "{}"'.format(
                str(time_tag_mjd_tandemx), str(time_tag_mjd_initial)
            )
        )
        logging.error(error_msg)
        raise ValueError(error_msg)

    # conversion step 2: fnf mask geotiff formatting: one single layer geotiff for each fnf "tamdem-x tile" (they will be merged togheter after equi7 conversion )
    fnf_mask_ground_dir_name = os.path.join(output_folder, 'fnf_mask_ground')
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
        equi7_fnf_mask_parent_tempdir = os.path.join(output_folder, 'fnf_mask_not_merged_tile_{}'.format(idx))
        equi7_out_name_curr = image2equi7grid(
            e7g,
            fnf_name_curr,
            equi7_fnf_mask_parent_tempdir,
            gdal_path=gdal_path,
            inband=None,
            subgrid_ids=None,
            accurate_boundary=False,
            resampling_type='bilinear',
            tile_nodata=np.float(0),
        )
        equi7_fnf_mask_tempdir_tiles_not_merged.extend(equi7_out_name_curr)

    # conversion step 4: merging fnf tiles
    equi7_fnf_mask_parent_tempdir = os.path.join(output_folder, 'initial_fnf_mask_equi7')
    ### managing final mask creation:
    equi7_fnf_mask_fnames = fnf_equi7_masks_merging(
        equi7_fnf_mask_tempdir_tiles_not_merged, equi7_fnf_mask_parent_tempdir
    )

    logging.info('...fnf mask conversion into equi7 done.\n')

    return equi7_fnf_mask_fnames


def fnf_equi7_load_filter_equi7format(
    forest_mask_catalogue_folder, e7g_curr_chain, product_resolution, output_folder, gdal_path
):

    # get the names of input equi7 fnf masks
    equi7_fnf_mask_tiff_names_list, equi7_subgrid_name = get_equi7_fnf_tiff_names(forest_mask_catalogue_folder)

    subgrid_code = equi7_subgrid_name[6:8]

    if not equi7_fnf_mask_tiff_names_list:
        error_str = 'Cannot find any Equi7 forest mask in input'
        logging.error(error_str)
        raise ValueError(error_str)

    # output folder of the mask that will be used in the current chainflow:
    equi7_fnf_mask_parent_dir = os.path.join(output_folder, 'initial_fnf_mask_equi7')

    # temporary folder for the input mask conversion
    equi7_fnf_mask_parent_tempdir = os.path.join(output_folder, 'temp_fnf_dir')

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
        intermediate_filtered_tiff_name = os.path.join(
            equi7_fnf_mask_parent_tempdir, 'initial_fnf_mask_equi7_{}.tif'.format(idx)
        )

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
                'Input FNF Equi7 mask geotransform is not in the correct format :'.format(subgrid_code) + str(e),
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
        error_str = 'Cannot find a valid input FNF Equi7 mask'
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
            temp_output_folder,
            'equi7_fnf_mask_tiles_merged',
            equi7_subgrid_name,
            equi7_tile_name,
            '_fnf_mask_' + equi7_subgrid_name + '_' + equi7_tile_name + '.tif',
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
                if not 'data_sum' in locals():
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
        outdata.FlushCache()  ##saves to disk!!
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
                temp_output_folder,
                'equi7_global_mask',
                stack_id,
                equi7_subgrid_name,
                equi7_tile_name,
                stack_id + '_final_mask_' + equi7_subgrid_name + '_' + equi7_tile_name + '.tif',
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
            if not 'out_mosaiked_data' in locals():
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
        out_mosaiked_equi7_name = os.path.join(
            out_tile_equi7_full_folder, 'FH_' + equi7_subgrid_folder_name + '_' + equi7_tile_name + '.tif'
        )

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
        outdata.FlushCache()  ##saves to disk!!
        outdata = None

        out_mosaiked_equi7_tiff_names.append(out_mosaiked_equi7_name)

        del out_mosaiked_data

    return out_mosaiked_equi7_tiff_names


def apply_calibration_screens(beta0_calibrated, beta0_raster_info, cal_screens, screens_raster_info, master_id):

    for idx, pf_name in enumerate(beta0_calibrated.keys()):

        if idx == 0:
            rg_axis_index = 0
            az_axis_index = 1

            # input original axis
            max_val_in_rg = screens_raster_info.num_samples * screens_raster_info.pixel_spacing_slant_rg
            rg_ax_in = np.arange(0, max_val_in_rg, raster_info.pixel_spacing_slant_rg)

            max_val_in_az = screens_raster_info.num_lines * screens_raster_info.pixel_spacing_az
            az_ax_in = np.arange(0, max_val_in_az, screens_raster_info.pixel_spacing_az)

            # output interpolated axis
            max_val_out_rg = beta0_raster_info.num_samples * beta0_raster_info.pixel_spacing_slant_rg
            rg_ax_out = np.arange(0, max_val_out_rg, beta0_raster_info.pixel_spacing_slant_rg)

            max_val_out_az = beta0_raster_info.num_lines * beta0_raster_info.pixel_spacing_az
            az_ax_out = np.arange(0, max_val_out_az, beta0_raster_info.pixel_spacing_az)

        # interpolate current screen over the data definiton axes:
        cal_scree_interp = cal_screens[pf_name]

        # interpolation of data along axis_index
        interp_fun = interp1d(rg_ax_in, cal_scree_interp, rg_axis_index, bounds_error=0)
        cal_scree_interp = interp_fun(rg_ax_out)
        interp_fun = interp1d(az_ax_in, cal_scree_interp, az_axis_index, bounds_error=0)
        cal_scree_interp = interp_fun(az_ax_out)

        for pol_id in beta0_calibrated[pf_name]:
            beta0_calibrated[pf_name][pol_id] = np.multiply(
                beta0_calibrated[pf_name][pol_id], np.exp(-1j * cal_scree_interp)
            )

    return beta0_calibrated


def apply_dem_flattening(beta0_calibrated, kz_in, reference_height, master_id, raster_info):

    for pf_name in beta0_calibrated.keys():
        for pol_id in beta0_calibrated[pf_name]:
            beta0_calibrated[pf_name][pol_id] = np.multiply(
                beta0_calibrated[pf_name][pol_id], np.exp(-1j * np.multiply(kz_in[pf_name], reference_height))
            )

    return beta0_calibrated


def format_folder_name():
    # date string format is: BIOMASS_L2_YYYYMMDDTHHMMSS
    now = datetime.now()
    year_str = str(now.year)
    month_str = str(now.month).rjust(2, '0')
    day_str = str(now.day).rjust(2, '0')
    h_str = str(now.hour).rjust(2, '0')
    m_str = str(now.minute).rjust(2, '0')
    s_str = str(now.second).rjust(2, '0')
    current_date_string = 'BIOMASS_L2_' + year_str + month_str + day_str + 'T' + h_str + m_str + s_str

    return current_date_string


def compute_and_oversample_geometry_auxiliaries(
    L1c_repository,
    aux_file_names,
    stack_id,
    acquisition_pf_names,
    master_id,
    enable_resampling,
    comp_ref_h=True,
    comp_dist=True,
    force_ellipsoid=False,
    sar_geometry_master=None,
):
    ecef_grid = None
    off_nadir_angle_rad = None
    kz = None
    reference_height = None
    R = None
    slope = None

    logging.info('Geometry library: loading auxiliary data: DEM...')
    pf_name = os.path.basename(aux_file_names.dem_folder)
    folder = os.path.dirname(aux_file_names.dem_folder)
    path_dem = os.path.join(folder, pf_name)

    try:
        pf_dem = ProductFolder(path_dem, 'r')
    except Exception as e:
        logging.error('Geometry library:  error during auxiliary DEM reading: ' + str(e), exc_info=True)
        raise

    ch_dem = pf_dem.get_channel(0)
    ri_dem = ch_dem.get_raster_info(0)

    dem_samples = ri_dem.samples
    lat_axis = (np.arange(0, ri_dem.lines) * ri_dem.lines_step + ri_dem.lines_start) * np.pi / 180
    lon_axis = (np.arange(0, dem_samples) * ri_dem.samples_step + ri_dem.samples_start) * np.pi / 180

    dem = pf_dem.read_data(0, [0, 0, ri_dem.lines, dem_samples])

    if force_ellipsoid:
        logging.info(
            'Geometry library: the ellipsoid altitude has been forced by input flag, only the slope over ellipsoid will be computed, after used to correcting geometry global / local reference.'
        )
        dem = 0 * dem

    logging.info('...done\n')

    logging.info('Geometry library: loading orbits from stack metadata...')
    ch_list = []
    dsi_list = []
    swath_id_list = []
    swath_info_list = []

    # read master alone to absure it will be first in the list
    pf_master = ProductFolder(os.path.join(L1c_repository, master_id), 'r')
    ch_master = pf_master.get_channel(0)
    ch_list.append(ch_master)
    ri_temp = ch_master.get_raster_info(0)
    if ri_temp.samples_step_unit == 'm':
        ri_master = convert_rasterinfo_meters_to_seconds(ri_temp)

        # ri_list.append( ri_master )
        num_samples = ri_temp.samples
        num_lines = ri_temp.lines
        lines_start = ri_temp.lines_start
        lines_start_unit = ri_temp.lines_start_unit
        lines_step = ri_temp.lines_step
        lines_step_unit = ri_temp.lines_step_unit
        samples_start = ri_temp.samples_start
        samples_start_unit = ri_temp.samples_start_unit
        samples_step = ri_temp.samples_step
        samples_step_unit = ri_temp.samples_step_unit

    else:
        ri_master = ri_temp

        ri_meters = convert_rasterinfo_seconds_to_meters(ri_temp)
        # ri_list.append( ri_master )
        num_samples = ri_meters.samples
        num_lines = ri_meters.lines
        lines_start = ri_meters.lines_start
        lines_start_unit = ri_meters.lines_start_unit
        lines_step = ri_meters.lines_step
        lines_step_unit = ri_meters.lines_step_unit
        samples_start = ri_meters.samples_start
        samples_start_unit = ri_meters.samples_start_unit
        samples_step = ri_meters.samples_step
        samples_step_unit = ri_meters.samples_step_unit

    sv_master = geometric_lib.cut_state_vectors(
        ch_master.get_state_vectors(0),
        ri_master.lines_start,
        ri_master.lines_start + ri_master.lines * ri_master.lines_step,
    )

    # Initiaize sar geometry
    flag_compute_sar_geometry = False
    if not sar_geometry_master:
        flag_compute_sar_geometry = True
        sar_geometry_master = geometric_lib.SARGeometry(ri_master, sv_master, lat_axis, lon_axis, dem)
        sar_geometry_master.compute_xyz()

    dsi_master = ch_master.get_dataset_info(0)
    dsi_list.append(dsi_master)

    swath_id_list.append(ch_master.get_swath_info(0).swath)
    swath_info_list.append(ch_master.get_swath_info(0))

    for acq_id in acquisition_pf_names:
        if not acq_id == master_id:
            pf_slave = ProductFolder(os.path.join(L1c_repository, acq_id), 'r')
            ch_slave = pf_slave.get_channel(0)
            ch_list.append(ch_slave)

            ri_temp = ch_slave.get_raster_info(0)
            if ri_temp.samples_step_unit == 'm':
                ri_slave = convert_rasterinfo_meters_to_seconds(ri_temp)
            else:
                ri_slave = ri_temp

            # ri_list. append( ri_slave )
            dsi_list.append(ch_slave.get_dataset_info(0))
            swath_id_list.append(ch_slave.get_swath_info(0).swath)
            swath_info_list.append(ch_slave.get_swath_info(0))

    logging.info('...done\n')

    logging.info('Geometry library: computing auxiliary data: slope...')

    if flag_compute_sar_geometry:
        sar_geometry_master.compute_terrain_slope()

    slope = sar_geometry_master.terrain_slope

    if enable_resampling:
        slope = data_oversample(slope, OVERSAMPLING_FACTOR, ri_master)[0]

    if force_ellipsoid:

        slope = slope.transpose()

        return ecef_grid, off_nadir_angle_rad, slope, kz, reference_height, R, sar_geometry_master

    else:
        pf_name = os.path.basename(aux_file_names.slope_file_names[stack_id])
        folder_name = os.path.dirname(aux_file_names.slope_file_names[stack_id])
        out_pf_name = os.path.join(folder_name, stack_id)

        if os.path.exists(out_pf_name):
            logging.warning('Overwriting of Auxiliary slope for stack ' + stack_id)
            shutil.rmtree(out_pf_name)
        if not os.path.exists(folder_name):
            os.makedirs(folder_name)

        # slope = sar_geometry_master.terrain_slope
        pf = ProductFolder(out_pf_name, 'w')

        # append a channel with the correct dimensions
        pf.append_channel(num_lines, num_samples, 'FLOAT32', header_offset=0)

        # write the ECEF coordinate to the channel
        pf.write_data(0, slope)

        # prepare the metadata elements
        data_channel_obj = pf.get_channel(0)
        metadata_obj = data_channel_obj.metadata
        metadatachannel_obj = metadata_obj.get_metadata_channels(0)

        # Raster Info
        ri = metadatachannel_obj.get_element('RasterInfo')
        ri.set_lines_axis(lines_start, lines_start_unit, lines_step, lines_step_unit)
        ri.set_samples_axis(samples_start, samples_start_unit, samples_step, samples_step_unit)

        # DataSetInfo: to insert the X, Y or Z description
        dsi_list[0].description = 'Auxiliary data: Slope [rad]'
        metadatachannel_obj.insert_element(dsi_list[0])

        pf.write_metadata(0)

        pf = None

    print('...done.\n')

    logging.info('Geometry library: computing auxiliary data: ECEFGRID...')

    ecef_grid = {
        'X': sar_geometry_master.x_sar_coords,
        'Y': sar_geometry_master.y_sar_coords,
        'Z': sar_geometry_master.z_sar_coords,
    }

    if enable_resampling:
        ecef_grid = data_oversample(ecef_grid, OVERSAMPLING_FACTOR, ri_master)[0]

    pf_name = os.path.basename(aux_file_names.ECEF_grid_file_names[stack_id])
    folder_name = os.path.dirname(aux_file_names.ECEF_grid_file_names[stack_id])
    out_pf_name = os.path.join(folder_name, stack_id)

    if os.path.exists(out_pf_name):
        logging.warning('Overwriting of Auxiliary ECEFGRID for stack ' + stack_id)
        shutil.rmtree(out_pf_name)
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)

    pf = ProductFolder(out_pf_name, 'w')
    dsi_ecefgrid = DataSetInfo('STRIPMAP', 435000000)
    dsi_ecefgrid.description = ''
    dsi_ecefgrid.image_type = 'SLC'
    dsi_ecefgrid.projection = 'SLANT RANGE'
    dsi_ecefgrid.side_looking = 'RIGHT'
    dsi_ecefgrid.sensor_name = 'BIOMASS'
    dsi_ecefgrid.sense_date = PreciseDateTime()
    dsi_ecefgrid.acquisition_station = 'NOT_AVAILABLE'
    dsi_ecefgrid.processing_center = 'NOT_AVAILABLE'
    dsi_ecefgrid.processing_date = PreciseDateTime()
    dsi_ecefgrid.processing_software = 'NOT_AVAILABLE'

    coord_names = ['X', 'Y', 'Z']
    for channel_idx, coord_name in enumerate(coord_names):

        pf.append_channel(num_lines, num_samples, 'FLOAT32', header_offset=0)
        pf.write_data(channel_idx, ecef_grid[coord_name])

        # prepare the metadata elements
        data_channel_obj = pf.get_channel(channel_idx)
        metadata_obj = data_channel_obj.metadata
        metadatachannel_obj = metadata_obj.get_metadata_channels(0)

        # Raster Info
        ri = metadatachannel_obj.get_element('RasterInfo')
        ri.set_lines_axis(lines_start, lines_start_unit, lines_step, lines_step_unit)
        ri.set_samples_axis(samples_start, samples_start_unit, samples_step, samples_step_unit)

        # DataSetInfo: to insert the X, Y or Z description
        dsi_ecefgrid.description = 'Auxiliary data: ' + coord_name + ' ECEF GRID [m]'
        metadatachannel_obj.insert_element(dsi_ecefgrid)

        pf.write_metadata(channel_idx)

    pf = None
    logging.info('...done\n')

    targets_coords = np.vstack(
        (
            sar_geometry_master.x_sar_coords.flatten(),
            sar_geometry_master.y_sar_coords.flatten(),
            sar_geometry_master.z_sar_coords.flatten(),
        )
    )

    if comp_ref_h:
        logging.info('Geometry library: computing auxiliary data: reference height...')
        reference_height = sar_geometry_master.dem_sar

        if enable_resampling:
            reference_height = data_oversample(reference_height, OVERSAMPLING_FACTOR, ri_master)[0]

        pf_name = os.path.basename(aux_file_names.reference_height_file_names[stack_id])
        folder_name = os.path.dirname(aux_file_names.reference_height_file_names[stack_id])
        out_pf_name = os.path.join(folder_name, stack_id)

        if os.path.exists(out_pf_name):
            logging.warning('Overwriting of Auxiliary reference height for stack ' + stack_id)
            shutil.rmtree(out_pf_name)
        if not os.path.exists(folder_name):
            os.makedirs(folder_name)

        # create the Product Folder in writing mode:
        pf = ProductFolder(out_pf_name, 'w')
        # append a channel with the correct dimensions
        pf.append_channel(num_lines, num_samples, 'FLOAT32', header_offset=0)
        # write the DTM to PF (reference height)
        pf.write_data(0, reference_height)

        # prepare the metadata elements
        data_channel_obj = pf.get_channel(0)
        metadata_obj = data_channel_obj.metadata
        metadatachannel_obj = metadata_obj.get_metadata_channels(0)

        # Raster Info
        ri = metadatachannel_obj.get_element('RasterInfo')
        ri.set_lines_axis(lines_start, lines_start_unit, lines_step, lines_step_unit)

        ri.set_samples_axis(samples_start, samples_start_unit, samples_step, samples_step_unit)

        # DataSetInfo
        dsi = DataSetInfo('STRIPMAP', 435000000)
        dsi.description = 'Auxiliary data: Reference height (DTM)'
        dsi.image_type = 'SLC'
        dsi.projection = 'SLANT RANGE'
        dsi.side_looking = 'RIGHT'
        dsi.sensor_name = 'BIOMASS'
        dsi.sense_date = PreciseDateTime()
        dsi.acquisition_station = 'NOT_AVAILABLE'
        dsi.processing_center = 'NOT_AVAILABLE'
        dsi.processing_date = PreciseDateTime()
        dsi.processing_software = 'NOT_AVAILABLE'

        metadatachannel_obj.insert_element(dsi)

        pf.write_metadata(0)
        pf = None

        print('...done.\n')

    if comp_dist:
        logging.info('Geometry library: computing auxiliary data: slant range distances...')

        pf_name = os.path.basename(aux_file_names.slant_range_distances_file_names[stack_id])
        folder_name = os.path.dirname(aux_file_names.slant_range_distances_file_names[stack_id])
        out_pf_name = os.path.join(folder_name, stack_id)

        if os.path.exists(out_pf_name):
            logging.warning('Overwriting of Auxiliary slant range distances for stack ' + stack_id)
            shutil.rmtree(out_pf_name)
        if not os.path.exists(folder_name):
            os.makedirs(folder_name)

        pf = ProductFolder(out_pf_name, 'w')
        R = {}
        for channel_idx, channel_curr in enumerate(ch_list):

            if channel_idx == 0:
                # distance between targets and master orbit
                R[swath_id_list[channel_idx]] = sar_geometry_master.compute_distance_master()
            else:

                # distance between targets and slave orbit
                ri_temp = channel_curr.get_raster_info(0)
                if ri_temp.samples_step_unit == 'm':
                    ri_slave = convert_rasterinfo_meters_to_seconds(ri_temp)
                else:
                    ri_slave = ri_temp

                sv_slave = geometric_lib.cut_state_vectors(
                    channel_curr.get_state_vectors(0),
                    ri_slave.lines_start,
                    ri_slave.lines_start + ri_slave.lines * ri_slave.lines_step,
                )

                gso_slave = create_general_sar_orbit(sv_slave)

                _, range_points_master_slave = geometric_lib.perform_inverse_geocoding(gso_slave, targets_coords)

                distance_slave = range_points_master_slave * LIGHTSPEED / 2
                distance_slave.shape = sar_geometry_master.x_sar_coords.shape

                R[swath_id_list[channel_idx]] = distance_slave

            if enable_resampling:
                R[swath_id_list[channel_idx]] = data_oversample(
                    R[swath_id_list[channel_idx]], OVERSAMPLING_FACTOR, ri_master
                )[0]

            # append a channel with the correct dimensions
            pf.append_channel(num_lines, num_samples, 'FLOAT32', header_offset=0)

            # write the R coordinate to the channel
            pf.write_data(channel_idx, R[swath_id_list[channel_idx]])

            # prepare the metadata elements
            data_channel_obj = pf.get_channel(channel_idx)
            metadata_obj = data_channel_obj.metadata
            metadatachannel_obj = metadata_obj.get_metadata_channels(0)

            # Raster Info
            ri = metadatachannel_obj.get_element('RasterInfo')
            ri.set_lines_axis(lines_start, lines_start_unit, lines_step, lines_step_unit)
            ri.set_samples_axis(samples_start, samples_start_unit, samples_step, samples_step_unit)

            # DataSetInfo: to insert the X, Y or Z description
            dsi_list[channel_idx].description = 'Auxiliary data: slant range distances [m]'
            metadatachannel_obj.insert_element(dsi_list[channel_idx])

            metadatachannel_obj.insert_element(swath_info_list[channel_idx])

            pf.write_metadata(channel_idx)
        pf = None

        print('...done.\n')

    logging.info('Geometry library: computing auxiliary data: off nadir...')

    pf_name = os.path.basename(aux_file_names.off_nadir_angle_file_names[stack_id])
    folder_name = os.path.dirname(aux_file_names.off_nadir_angle_file_names[stack_id])
    out_pf_name = os.path.join(folder_name, stack_id)

    if os.path.exists(out_pf_name):
        logging.warning('Overwriting of Auxiliary off nadir for stack ' + stack_id)
        shutil.rmtree(out_pf_name)
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)

    pf = ProductFolder(out_pf_name, 'w')
    off_nadir_angle_rad = {}
    for channel_idx, channel_curr in enumerate(ch_list):

        if channel_idx == 0:
            # master
            sensor_positions_master = sar_geometry_master.get_sensor_positions_master()
            off_nadir_angle_rad[
                swath_id_list[channel_idx]
            ] = sar_geometry_master.compute_off_nadir_angles_from_positions(
                sensor_positions_master, sensor_positions_master
            )
        else:
            # slaves
            ri_temp = channel_curr.get_raster_info(0)
            if ri_temp.samples_step_unit == 'm':
                ri_slave = convert_rasterinfo_meters_to_seconds(ri_temp)
            else:
                ri_slave = ri_temp

            sv_slave = geometric_lib.cut_state_vectors(
                channel_curr.get_state_vectors(0),
                ri_slave.lines_start,
                ri_slave.lines_start + ri_slave.lines * ri_slave.lines_step,
            )

            gso_slave = create_general_sar_orbit(sv_slave)

            azimuth_points_master_slave, _ = geometric_lib.perform_inverse_geocoding(gso_slave, targets_coords)

            sensor_positions_slave = gso_slave.get_position(azimuth_points_master_slave)
            sensor_positions_slave.shape = (3, *sar_geometry_master.x_sar_coords.shape)
            off_nadir_angle_rad[
                swath_id_list[channel_idx]
            ] = sar_geometry_master.compute_off_nadir_angles_from_positions(
                sensor_positions_master, sensor_positions_slave
            )

        if enable_resampling:
            off_nadir_angle_rad[swath_id_list[channel_idx]] = data_oversample(
                off_nadir_angle_rad[swath_id_list[channel_idx]], OVERSAMPLING_FACTOR, ri_master
            )[0]

        # append a channel with the correct dimensions
        pf.append_channel(num_lines, num_samples, 'FLOAT32', header_offset=0)

        # write the ECEF coordinate to the channel
        pf.write_data(channel_idx, off_nadir_angle_rad[swath_id_list[channel_idx]])

        # prepare the metadata elements
        data_channel_obj = pf.get_channel(channel_idx)
        metadata_obj = data_channel_obj.metadata
        metadatachannel_obj = metadata_obj.get_metadata_channels(0)

        # Raster Info
        ri = metadatachannel_obj.get_element('RasterInfo')
        ri.set_lines_axis(lines_start, lines_start_unit, lines_step, lines_step_unit)
        ri.set_samples_axis(samples_start, samples_start_unit, samples_step, samples_step_unit)

        # DataSetInfo: to insert the X, Y or Z description
        dsi_list[channel_idx].description = 'Auxiliary data: Incidence Angle [rad]'
        metadatachannel_obj.insert_element(dsi_list[channel_idx])

        metadatachannel_obj.insert_element(swath_info_list[channel_idx])

        pf.write_metadata(channel_idx)

    pf = None
    print('...done.\n)')

    logging.info('Geometry library: computing auxiliary data: KZ (wavenumbers)...')

    pf_name = os.path.basename(aux_file_names.kz_file_names[stack_id])
    folder_name = os.path.dirname(aux_file_names.kz_file_names[stack_id])
    out_pf_name = os.path.join(folder_name, stack_id)

    if os.path.exists(out_pf_name):
        logging.warning('Overwriting of Auxiliary KZ (wavenumbers) for stack ' + stack_id)
        shutil.rmtree(out_pf_name)
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)

    pf = ProductFolder(out_pf_name, 'w')
    kz = {}
    for channel_idx, channel_curr in enumerate(ch_list):

        if channel_idx == 0:
            # Master
            wavelength = LIGHTSPEED / ch_master.get_dataset_info(0).fc_hz

            kz[swath_id_list[channel_idx]] = np.zeros((ri_master.lines, ri_master.samples))
        else:
            kz[swath_id_list[channel_idx]] = geometric_lib.compute_vertical_wavenumber_angles(
                wavelength, off_nadir_angle_rad[swath_id_list[0]], off_nadir_angle_rad[swath_id_list[channel_idx]]
            )

        if enable_resampling:
            kz[swath_id_list[channel_idx]] = data_oversample(
                kz[swath_id_list[channel_idx]], OVERSAMPLING_FACTOR, ri_master
            )[0]

        # append a channel with the correct dimensions
        pf.append_channel(num_lines, num_samples, 'FLOAT32', header_offset=0)

        # write the KZ  to the channel
        pf.write_data(channel_idx, kz[swath_id_list[channel_idx]])

        # prepare the metadata elements
        data_channel_obj = pf.get_channel(channel_idx)
        metadata_obj = data_channel_obj.metadata
        metadatachannel_obj = metadata_obj.get_metadata_channels(0)

        # Raster Info
        ri = metadatachannel_obj.get_element('RasterInfo')
        ri.set_lines_axis(lines_start, lines_start_unit, lines_step, lines_step_unit)
        ri.set_samples_axis(samples_start, samples_start_unit, samples_step, samples_step_unit)

        # DataSetInfo: to insert the X, Y or Z description
        dsi_list[channel_idx].description = 'Auxiliary data: KZ (Interferometric wavenumber matrix)'
        metadatachannel_obj.insert_element(dsi_list[channel_idx])

        metadatachannel_obj.insert_element(swath_info_list[channel_idx])

        pf.write_metadata(channel_idx)

    pf = None
    logging.info('...done\n')

    # Transpose all the outputs as done in the case of auxiliaries directly read (and not computed), see "read_and_oversample_aux_data"
    if not ecef_grid is None:
        for coord_name in coord_names:
            ecef_grid[coord_name] = ecef_grid[coord_name].transpose()

    for swath_id in swath_id_list:
        if not off_nadir_angle_rad is None:
            off_nadir_angle_rad[swath_id] = off_nadir_angle_rad[swath_id].transpose()
        if not kz is None:
            kz[swath_id] = kz[swath_id].transpose()
        if not R is None:
            R[swath_id] = R[swath_id].transpose()

    if not slope is None:
        slope = slope.transpose()

    if not reference_height is None:
        reference_height = reference_height.transpose()

    return ecef_grid, off_nadir_angle_rad, slope, kz, reference_height, R, sar_geometry_master


def check_fnf_folder_format(fnf_mask_catalogue_folder):

    fnf_content_folders = os.listdir(fnf_mask_catalogue_folder)
    if not fnf_content_folders:
        error_msg = 'input fnf mask is mandatory'
        logging.error(error_msg)
        raise Exception(error_msg)

    found_tandemx = False
    found_equi7 = False
    for fnf_content_folder in fnf_content_folders:
        if 'TDM_FNF_' in fnf_content_folder:
            found_tandemx = True

        elif 'EQUI7_' in fnf_content_folder:
            found_equi7 = True

        else:
            error_msg = 'Cannot recognize input FNF mask format: valid formats are EQUI7 (an "EQUI7_*/" folder containing equi7 tiles sub-folders, example: "EQUI7_AF025M/E042N048T3) or TANDEM-X (collection of "TDM_FNF_*/" folders '
            logging.error(error_msg)
            raise Exception(error_msg)

    if found_tandemx and found_equi7:
        error_msg = ' Found bowth TANDEM-X and EQUI7 format FNF masks, this is not supported.'
        logging.error(error_msg)
        raise Exception(error_msg)
    elif found_tandemx:
        return 'TANDEM-X'
    elif found_equi7:
        return 'EQUI7'


def check_cal_format(reference_agb_folder):

    cal_format = 'RASTER'

    return cal_format


def get_equi7_fnf_tiff_names(fnf_mask_catalogue_folder):

    equi7_tiles_subfolders = []
    e7g_type = ''
    for idx, equi7_subgrid_name in enumerate(os.listdir(fnf_mask_catalogue_folder)):

        if idx == 1:
            logging.warning(
                'Found more than one sub-grid in input Equi7 format FNF mask: prototype supports just one grid at a time, only {} sub-grid will be read'.format(
                    equi7_subgrid_name
                )
            )
            break

        equi7_folder = os.path.join(fnf_mask_catalogue_folder, equi7_subgrid_name)
        temp_list = os.listdir(equi7_folder)
        if e7g_type == '':
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

        cal_names.append(os.path.join(reference_agb_folder, tiff_name))

    return cal_names


def get_foss_cal_names(reference_agb_folder):

    tiff_names = os.listdir(reference_agb_folder)

    cal_names = []
    for tiff_name in tiff_names:

        cal_names.append(os.path.join(reference_agb_folder, tiff_name))

    return cal_names


def get_min_time_stamp_repository(L1c_repository, stack_composition):

    for idx, (unique_stack_id, unique_acq_pf_names) in enumerate(stack_composition.items()):

        if idx == 0:
            time_tag_mjd_initial = get_min_time_stamp_repository_core(L1c_repository, unique_acq_pf_names)
        else:
            time_tag_mjd_initial = min(
                time_tag_mjd_initial, get_min_time_stamp_repository_core(L1c_repository, unique_acq_pf_names)
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


def collect_stacks_to_be_merged(stack_composition):

    stacks_to_merge_dict = {}
    for unique_stack_id in stack_composition.keys():

        (
            global_cycle_idx,
            heading_deg,
            rg_swath_idx,
            rg_sub_swath_idx,
            az_swath_idx,
            baseline_idx,
        ) = decode_unique_acquisition_id_string(stack_composition[unique_stack_id][0], output_format='string')

        unique_merged_stack_id = (
            'GC_' + global_cycle_idx + '_RGSW_' + rg_swath_idx + '_RGSBSW_' + rg_sub_swath_idx + '_AZSW_' + az_swath_idx
        )

        if unique_merged_stack_id in stacks_to_merge_dict.keys():
            stacks_to_merge_dict[unique_merged_stack_id].append(unique_stack_id)
        else:
            stacks_to_merge_dict[unique_merged_stack_id] = [unique_stack_id]

    return stacks_to_merge_dict


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
        error_str = 'Stack heading not valid: {}  [deg]'.format(heading_deg)
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


def convert_rasterinfo_meters_to_seconds(ri_meters):

    ri_seconds = RasterInfo(
        ri_meters.lines,
        ri_meters.samples,
        ri_meters.cell_type,
        ri_meters.file_name,
        ri_meters.header_offset_bytes,
        ri_meters.row_prefix_bytes,
        ri_meters.byte_order,
    )

    lines_step_s = ri_meters.lines_step / SATELLITE_VELOCITY
    lines_step_unit = 's'
    samples_start_s = ri_meters.samples_start / LIGHTSPEED * 2
    samples_start_unit = 's'
    samples_step_s = ri_meters.samples_step / LIGHTSPEED * 2
    samples_step_unit = 's'

    ri_seconds.set_lines_axis(ri_meters.lines_start, ri_meters.lines_start_unit, lines_step_s, lines_step_unit)
    ri_seconds.set_samples_axis(samples_start_s, samples_start_unit, samples_step_s, samples_step_unit)

    return ri_seconds


def convert_rasterinfo_seconds_to_meters(ri_seconds):

    ri_meters = RasterInfo(
        ri_seconds.lines,
        ri_seconds.samples,
        ri_seconds.cell_type,
        ri_seconds.file_name,
        ri_seconds.header_offset_bytes,
        ri_seconds.row_prefix_bytes,
        ri_seconds.byte_order,
    )

    lines_step_m = ri_meters.lines_step * SATELLITE_VELOCITY
    lines_step_unit = 'm'
    samples_start_m = ri_meters.samples_start * LIGHTSPEED / 2
    samples_start_unit = 'm'
    samples_step_m = ri_meters.samples_step * LIGHTSPEED / 2
    samples_step_unit = 'm'

    ri_meters.set_lines_axis(ri_seconds.lines_start, ri_seconds.lines_start_unit, lines_step_m, lines_step_unit)
    ri_meters.set_samples_axis(samples_start_m, samples_start_unit, samples_step_m, samples_step_unit)

    return ri_meters
