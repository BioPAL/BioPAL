# SPDX-FileCopyrightText: BioPAL <biopal@esa.int>
# SPDX-License-Identifier: MIT

import os
import logging
import collections
from shapely.geometry import MultiPoint
from biopal.utility.utility_functions import getBinaryNameFromChannelIDX
from biopal.io.xml_io import (
    parse_input_file,
    geographic_boundaries_latlon,
    write_input_file,
    stack_based_processing,
)
from biopal.utility.utility_functions import Task, start_logging
from biopal.io.data_io import readBiomassHeader_core
from arepytools.timing.precisedatetime import PreciseDateTime


class InvalidInputError(ValueError):
    pass


class dataset_query(Task):
    """
    Performs a query on the BioPAL dataSet: the query rule is specified 
    into the input file "dataset_query" section to get stacks from a 
    specified geographical and temporal region.
    Writes in output an updated input file, filled with new section 
    "stack_based_processing"
    
    Parameters
    ----------
    input_file : BioPAL xml input file path, or input file object
        The input file shoud contain mandatory following sections:
        L2_product, output_specification, dataset_query.
        See as an example: BioPAL\inputs\Input_File.xml.
        input_file can alternatively be an object, parsed from an xml file as:
            from biopal.io.xml_io import parse_input_file
            input_file = parse_input_file(input_file_xml)
        
    Returns
    -------
    input_file_up_xml
        BioPAL xml input file path, updated with "stack_based_processing" new section
    stack_composition
        dictionary containing, for each stack ID key, its acquisition names;
        info contained also on the stack_based_processing new section
    geographic_boundaries
        Minimum and maximum latitudes and longitudes box which contains all the queried stacks;
        info contained also on the stack_based_processing new section
        
    Examples
    -------
    >>> input_path = BioPAL / "inputs" / "Input_File.xml"
    >>> from biopal.agb.main_AGB import dataset_query
    >>> dataset_query_obj = dataset_query()
    >>> input_file_up_xml = dataset_query.run( input_path )
    >>> print('The Input_File has been updated with the new section "stack_based_processing" and saved to {}'.format(input_file_up_xml))
    """

    def __init__(self):
        super().__init__()

    def _run(self, input_file):

        if isinstance(input_file, str):
            input_params_obj = parse_input_file(input_file)
        else:
            input_params_obj = input_file
        if (
            input_params_obj.L2_product is None
            or input_params_obj.output_specification is None
            or input_params_obj.dataset_query is None
        ):

            error_message = [
                "BioPAL APP input file should contain at least"
                ' "L2_product", "output_specification" and "dataset_query" sections'
            ]
            raise InvalidInputError(error_message)

        if not hasattr(logging.getLoggerClass().root.handlers[0], "baseFilename"):
            start_logging(
                input_params_obj.output_specification.output_folder,
                input_params_obj.L2_product,
                "DEBUG",
                app_name="dataset_query",
            )
        (
            stack_composition,
            geographic_boundaries,
            geographic_boundaries_per_stack,
        ) = data_select_by_date_and_boundaries(input_params_obj.dataset_query)
        if not len(stack_composition):
            logging.error(
                "Cannot find any data which meets the conditions to run the processing: see the log warnings. \n"
            )
            logging.info("BIOMASS L2 Processor ended without any processing performed.")
            return 0

        logging.info("Geographic Boundaries extracted [deg]:")
        logging.info("Latitude  min: " + str(geographic_boundaries.lat_min))
        logging.info("Latitude  max: " + str(geographic_boundaries.lat_max))
        logging.info("Longitude min: " + str(geographic_boundaries.lon_min))
        logging.info("Longitude max: " + str(geographic_boundaries.lon_max))

        # write the input files for the activated chains:
        stack_based_proc_obj = fill_stack_based_processing_obj(
            input_params_obj, stack_composition, geographic_boundaries, geographic_boundaries_per_stack,
        )

        input_file_up_xml = os.path.join(
            input_params_obj.output_specification.output_folder,
            "InputFile_StackBasedProcessing{}.xml".format(input_params_obj.L2_product),
        )
        input_params_obj.output_specification.output_folder = os.path.join(
            input_params_obj.output_specification.output_folder, input_params_obj.L2_product,
        )
        input_params_obj.stack_based_processing = stack_based_proc_obj
        write_input_file(input_params_obj, input_file_up_xml)

        end_message = ("Query completed.")
        logging.info(end_message)
        print(end_message)
        
        return input_file_up_xml


def data_select_by_date_and_boundaries(dataset_query_obj):

    geographic_boundaries_curr = geographic_boundaries_latlon(None, None, None, None)

    # from all the database of data, selectonly the ones which matches the
    # input parameters and return a dictionary with the stacks to be used

    start_date = dataset_query_obj.L1C_date.start
    stop_date = dataset_query_obj.L1C_date.stop

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
    pf_names_list = os.listdir(dataset_query_obj.L1C_repository)

    for pf_idx, pf_name in enumerate(pf_names_list):
        date_is_matching = False
        boundaries_are_matching = False
        stack_is_matching = False
        num_acq_are_enough = False

        if len(pf_name) != 47:
            error_msg = 'Unrecognized uniqie acquisition ID (it will be skipped): "{}" from repository {}; the correct format is "GC_XX_H_XXX.XX_RGSW_XX_RGSBSW_X_AZSW_XX_BSL_XX" where each "X" is an integer'.format(
                pf_name, dataset_query_obj.L1C_repository
            )
            logging.error(error_msg)
            continue

        output_pf_full_name = os.path.join(dataset_query_obj.L1C_repository, pf_name)

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
        (
            datestr,
            lon_min,
            lon_max,
            lat_min,
            lat_max,
            unique_acq_id,
            resolution_m_slant_rg,
            resolution_m_az,
            sensor_velocity,
        ) = readBiomassHeader_core(raster_file)

        uniqie_stack_id = unique_acq_id[0:40]

        # first test:check the date:
        date = PreciseDateTime().set_from_utc_string(datestr)

        if pf_idx == 0:
            edge_values = {
                "date_min": date,
                "date_max": date,
                "lon_min": lon_min,
                "lon_max": lon_max,
                "lat_min": lat_min,
                "lat_max": lat_max,
            }

        else:
            edge_values["date_min"] = min(edge_values["date_min"], date)
            edge_values["date_max"] = max(edge_values["date_max"], date)
            edge_values["lon_min"] = min(edge_values["lon_min"], lon_min)
            edge_values["lon_max"] = max(edge_values["lon_max"], lon_max)
            edge_values["lat_min"] = min(edge_values["lat_min"], lat_min)
            edge_values["lat_max"] = max(edge_values["lat_max"], lat_max)

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
            geographic_boundaries_per_stack = geographic_boundaries_latlon(None, None, None, None)
            geographic_boundaries_per_stack.lon_min = lon_min
            geographic_boundaries_per_stack.lon_max = lon_max
            geographic_boundaries_per_stack.lat_min = lat_min
            geographic_boundaries_per_stack.lat_max = lat_max
            stacks_geographic_boundaries[uniqie_stack_id] = geographic_boundaries_per_stack

        # check 2/4 search for data inside the user dates boundaries:
        date_is_matching = check_on_dates(start_date, stop_date, date)

        # check 3/4 search for data inside the user ground coordinates boundaries:
        boundaries_are_matching = check_on_boundaries(
            lon_min, lon_max, lat_min, lat_max, dataset_query_obj.boundaries.point
        )

        # check 4/4 (optional) ckeck if data matches one of the user specified stack:
        stack_is_matching = check_on_stacks_to_find(uniqie_stack_id, dataset_query_obj.stacks_to_find)

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
        logging.warning("\n")
        logging.warning("Cannot find any valid data:")
        if len(stacks_composition_only_dates_match) or len(stacks_composition_only_boundaries_match):

            logging.warning("    be aware that data In L1cRepository are confined inside following boundaries:")
            logging.warning("    minimun date is       {}".format(edge_values["date_min"]))
            logging.warning("    maximun date is       {}".format(edge_values["date_max"]))
            logging.warning("    minimun Longitude is  {}".format(edge_values["lon_min"]))
            logging.warning("    maximum Longitude is  {}".format(edge_values["lon_max"]))
            logging.warning("    minimun Latitude  is  {}".format(edge_values["lat_min"]))
            logging.warning("    maximum Latitude  is  {} \n".format(edge_values["lat_max"]))
            logging.warning("    Instead the user asked for following boundaries:")
            logging.warning("    user min date         {}".format(start_date))
            logging.warning("    user max date         {}".format(stop_date))

            point_dict = dataset_query_obj.boundaries.point

            logging.warning("    user min Longitude    {}".format(min([x["Longitude"] for x in point_dict])))
            logging.warning("    user max Longitude    {}".format(max([x["Longitude"] for x in point_dict])))
            logging.warning("    user min Latitude     {}".format(min([x["Latitude"] for x in point_dict])))
            logging.warning("    user max Latitude     {} \n".format(max([x["Latitude"] for x in point_dict])))

            if len(stacks_composition_only_dates_match):
                logging.warning(
                    "    Following stacks are confined in the desired L1cDates but not in the desired geographic boundaries"
                )
                logging.warning("    {} \n".format(list(stacks_composition_only_dates_match.keys())))

            if len(stacks_composition_only_boundaries_match):
                logging.warning(
                    "    Following stacks are confined in the desired geographic boundaries but not in the desired L1cDates"
                )
                logging.warning("    {} \n".format(list(stacks_composition_only_boundaries_match.keys())))

        if len(stacks_composition_not_enough_acq):
            logging.warning(
                "Following stacks contain just ONE acquisition (a minimum of TWO are needed to have ONE baseline)"
            )
            logging.warning("{} \n".format(list(stacks_composition_not_enough_acq.keys())))

        if (
            not len(stacks_composition_only_dates_match)
            and not len(stacks_composition_only_boundaries_match)
            and not len(stacks_composition_not_enough_acq)
        ):
            logging.warning("Cannot find any data, the L1cRepository is empty")
            pf_idx = 0

    else:
        logging.info("Found #{} stacks: {}".format(num_stacks, list(stacks_composition.keys())))
        logging.info("Stacks composition: \n")
        for stack_curr, pf_list in stacks_composition.items():
            pf_list.sort()
            logging.info("Stack " + stack_curr + " composition ( unique acquisitions id ):")
            for pf_name in pf_list:
                logging.info("    " + pf_name)

    stacks_composition = collections.OrderedDict(sorted(stacks_composition.items()))

    return stacks_composition, geographic_boundaries_curr, stacks_geographic_boundaries


def fill_stack_based_processing_obj(
    input_params_obj, stack_composition, geographic_boundaries, geographic_boundaries_per_stack,
):

    """
    Internal function called by the dataset_query APP:
        
        after the query over the dataSet, the datset_query APP fills the structure
        to be written into the xml input file for the next APP, which is the 
        stack based processing.
        The output structure contains all the queried stacks paths and names from dataSet ( 
        and auxiliary paths), other than boundaries.
        The returned object "stack_based_processing_obj" structure will have the specific 
        ields needed by the activated chain 
        (AGB, FD, FH or TOMO_FH), other than all the common fields to each chain.

        Usage of the returned object:
            The returned object can be added to the input_params_obj and it 
            can be written to disk if needed, as (stack_based_processing_obj is 
                                                  overwritten by this command, if already present):
                - input_params_obj.stack_based_processing_obj = stack_based_processing_obj
                - write_input_file(input_params_obj, input_file_xml)
            
    """

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
    ) = get_auxiliary_paths(stack_composition, input_params_obj.dataset_query.L1C_aux_repository)

    if input_params_obj.L2_product == "AGB":
        system_decorrelation_fun_folder = None
        average_covariance_folder = None

    if input_params_obj.L2_product == "FD":
        reference_agb_folder = None
        system_decorrelation_fun_folder = None
        forest_height_folder = None

    if input_params_obj.L2_product == "FH":
        reference_agb_folder = None
        average_covariance_folder = None
        forest_height_folder = None

    if input_params_obj.L2_product == "TOMO_FH":

        stack_composition_TOMO = {}
        for stack_id, acquisitions_list in stack_composition.items():
            if len(acquisitions_list) >= 3:
                stack_composition_TOMO[stack_id] = stack_composition[stack_id]

        len_stacks_all = stack_composition.keys()
        len_stacks_tomo = stack_composition_TOMO.keys()
        if len_stacks_tomo:
            if len_stacks_tomo < len_stacks_all:
                logging.info(
                    "The stacks with less than #3 acquisitions cannot be used to perform TOMO FH estimation \
                     and they will be removed: #{} stacks among #{} totaL, have been removed".format(
                        len_stacks_tomo, len_stacks_all
                    )
                )
        else:
            error_message = (
                "cannot find any stack with more than #2 acquisitions: TOMO FH estimation cannot be executed"
            )
            logging.error(error_message)
            raise InvalidInputError(error_message)

        stack_composition = stack_composition_TOMO
        reference_agb_folder = None
        forest_mask_catalogue_folder = None
        system_decorrelation_fun_folder = None
        average_covariance_folder = None
        forest_height_folder = None

    stack_based_processing_obj = stack_based_processing(
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
        geographic_boundaries,
        geographic_boundaries_per_stack,
    )

    return stack_based_processing_obj


def get_auxiliary_paths(stack_composition, L1c_aux_data_repository):
    # get the names of the folders containing of the auxiliary data: the names
    # shouls match the default ones

    ECEFGRID_folder = os.path.join(L1c_aux_data_repository, "Geometry", "ECEFGRID")
    KZ_folder = os.path.join(L1c_aux_data_repository, "Geometry", "KZ")
    off_nadir_angle_folder = os.path.join(L1c_aux_data_repository, "Geometry", "OffNadirAngles")
    reference_height_folder = os.path.join(L1c_aux_data_repository, "Geometry", "ReferenceHeight")
    slant_range_distances = os.path.join(L1c_aux_data_repository, "Geometry", "SlantRangeDistances")
    slope_folder = os.path.join(L1c_aux_data_repository, "Geometry", "Slope")

    forest_mask_catalogue_folder = os.path.join(L1c_aux_data_repository, "ForestMask")
    dem_folder = os.path.join(L1c_aux_data_repository, "DEM")
    reference_agb_folder = os.path.join(L1c_aux_data_repository, "ReferenceAGB")
    average_covariance_folder = os.path.join(L1c_aux_data_repository, "AverageCovariance")
    system_decorrelation_fun_folder = os.path.join(L1c_aux_data_repository, "SystemDecorrelationFunction")
    calibration_screens_folder = os.path.join(L1c_aux_data_repository, "CalibrationScreens")
    forest_height_folder = os.path.join(L1c_aux_data_repository, "ForestHeight")

    # Those are stack based auxiliary:
    ECEF_grid_file_names = {}
    kz_file_names = {}
    off_nadir_angle_file_names = {}
    reference_height_file_names = {}
    slant_range_distances_file_names = {}
    slope_file_names = {}
    calibration_screens_file_names = {}

    for stack_id in stack_composition.keys():

        ECEF_grid_file_names[stack_id] = os.path.join(r"" + ECEFGRID_folder, stack_id)
        kz_file_names[stack_id] = os.path.join(r"" + KZ_folder, stack_id)
        off_nadir_angle_file_names[stack_id] = os.path.join(r"" + off_nadir_angle_folder, stack_id)
        reference_height_file_names[stack_id] = os.path.join(r"" + reference_height_folder, stack_id)
        slant_range_distances_file_names[stack_id] = os.path.join(r"" + slant_range_distances, stack_id)
        slope_file_names[stack_id] = os.path.join(r"" + slope_folder, stack_id)
        calibration_screens_file_names[stack_id] = os.path.join(r"" + calibration_screens_folder, stack_id)

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
    data_longitudes = [x["Longitude"] for x in point_dict]
    data_latitudes = [x["Latitude"] for x in point_dict]
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
