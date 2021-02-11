"""biopal processor

Usage:
  biopal [--conf=CONFFOLDER] INPUTFILEXML
  
Arguments:
  INPUTFILEXML       Path of the xml input file
  
Options:
  --conf=CONFFOLDER  Set the directory containing the xml configuration files
  -h --help          Show this screen
  --version          Show version
"""

import os
import pkg_resources
import logging
import collections
from shutil import copyfile, which
from shapely.geometry import MultiPoint
from biopal.utility.utility_functions import (
    check_if_path_exists,
    format_folder_name,
    getBinaryNameFromChannelIDX,
    decode_unique_acquisition_id_string,
    set_gdal_paths,
)
from biopal.io.xml_io import (
    parse_biomassL2_main_input_file,
    parse_biopal_configuration_file,
    geographic_boundaries,
    write_chains_input_file,
    input_params,
    write_boundaries_files,
)
from biopal.io.data_io import readBiomassHeader_core

from arepytools.timing.precisedatetime import PreciseDateTime

# following imports are tested with "try" because different modalities of
# software releases may lack of some libraries, for example when chains are
# delivered separately
try:
    from biopal.agb.main_AGB import AboveGroundBiomass
except:
    pass
try:
    from biopal.fh.main_FH import ForestHeight
except:
    pass
try:
    from biopal.tomo_fh.main_TOMO_FH import TomoForestHeight
except:
    pass
try:
    from biopal.fd.main_FD import ForestDisturbance
except:
    pass
try:
    from biopal.tomo.main_TOMO_CUBE import main_TOMO_CUBE
except:
    pass


# main biopal:
# 1) Sets enviriment and logging, parses the main input
# 2) Chooses data to be processed
# 3) Launches each activated chain
def biomassL2_processor_run(input_file_xml, conf_folder=None):

    if conf_folder is None:
        biopal_configuration_file_xml = pkg_resources.resource_filename("biopal", "conf/biopal_configuration.xml")
        configuration_file_AGB = pkg_resources.resource_filename("biopal", "conf/ConfigurationFile_AGB.xml")
        configuration_file_FH = pkg_resources.resource_filename("biopal", "conf/ConfigurationFile_FH.xml")
        configuration_file_TOMO_FH = pkg_resources.resource_filename("biopal", "conf/ConfigurationFile_TOMO_FH.xml")
        configuration_file_FD = pkg_resources.resource_filename("biopal", "conf/ConfigurationFile_FD.xml")
    else:
        default_configuration_folder = conf_folder
        check_if_path_exists(default_configuration_folder, "FOLDER")

        biopal_configuration_file_xml = os.path.join(default_configuration_folder, "biopal_configuration.xml")
        configuration_file_AGB = os.path.join(default_configuration_folder, "ConfigurationFile_AGB.xml")
        configuration_file_FH = os.path.join(default_configuration_folder, "ConfigurationFile_FH.xml")
        configuration_file_TOMO_FH = os.path.join(default_configuration_folder, "ConfigurationFile_TOMO_FH.xml")
        configuration_file_FD = os.path.join(default_configuration_folder, "ConfigurationFile_FD.xml")
    
    # read the main input file
    main_input_struct = parse_biomassL2_main_input_file(input_file_xml)

    # date string format is: BIOMASS_L2_YYYYMMDDTHHMMSS
    current_date_string = format_folder_name()

    output_folder = os.path.join(main_input_struct.output_folder, current_date_string)

    if os.path.exists(output_folder):
        error_message = "Output Folder " + output_folder + " already exists, exiting now."
        raise RuntimeError(error_message)
    os.makedirs(output_folder)

    # start the main logging
    log_file_name = start_logging(output_folder, main_input_struct.proc_flags, "DEBUG")

    gdal_path, gdal_environment_path = parse_biopal_configuration_file(biopal_configuration_file_xml)
    gdal_path, gdal_environment_path = set_gdal_paths(gdal_path, gdal_environment_path)

    logging.debug("Configuration   folder is {}".format(conf_folder))
    logging.info(" Auxiliary data folder is {}".format(main_input_struct.L1c_aux_data_repository))
    logging.info("Results will be saved into output folder {}".format(output_folder))

    # get all the data that matches input datations and ground boundaries
    # a dictionary with stack_id and scene oroduct folder names is retrived
    logging.info("Searching data from data set : " + main_input_struct.L1c_repository)
    logging.info("Research will be done according to user dates and boundaries....\n")
    try:

        (
            stack_composition,
            geographic_boundaries,
            geographic_boundaries_per_stack,
        ) = data_select_by_date_and_boundaries(main_input_struct)
        if not len(stack_composition):
            logging.error(
                "Cannot find any data which meets the conditions to run the processing: see the log warnings. \n"
            )
            logging.info("BIOMASS L2 Processor ended without any processing performed.")
            return 0

    except Exception as e:
        logging.error(e, exc_info=True)
        raise

    write_boundaries_files(geographic_boundaries, geographic_boundaries_per_stack, output_folder)

    logging.info("Geographic Boundaries extracted [deg]:")
    logging.info("Latitude  min: " + str(geographic_boundaries.lat_min))
    logging.info("Latitude  max: " + str(geographic_boundaries.lat_max))
    logging.info("Longitude min: " + str(geographic_boundaries.lon_min))
    logging.info("Longitude max: " + str(geographic_boundaries.lon_max))

    # write the input files for the activated chains:
    (
        AGB_input_file_xml,
        FD_input_file_xml,
        FH_input_file_xml,
        TOMO_FH_input_file_xml,
        TOMO_input_file_xml,
    ) = write_chains_input_file_main(output_folder, main_input_struct, stack_composition)

    # Execute all the activated chains, separately

    if main_input_struct.proc_flags.FH or main_input_struct.proc_flags.TOMO_FH:
        stacks_to_merge_dict = collect_stacks_to_be_merged(stack_composition)

    # AGB
    if main_input_struct.proc_flags.AGB:

        agb_obj = AboveGroundBiomass(
            configuration_file_AGB, geographic_boundaries, geographic_boundaries_per_stack, gdal_path,
        )

        agb_obj.run(AGB_input_file_xml)

    # FH
    if main_input_struct.proc_flags.FH:

        fh_obj = ForestHeight(configuration_file_FH, geographic_boundaries, stacks_to_merge_dict, gdal_path,)

        fh_obj.run(FH_input_file_xml)

    # TOMO FH
    if main_input_struct.proc_flags.TOMO_FH:

        tomo_fh_obj = TomoForestHeight(configuration_file_TOMO_FH, stacks_to_merge_dict, gdal_path,)

        tomo_fh_obj.run(TOMO_FH_input_file_xml)

    # FD
    if main_input_struct.proc_flags.FD:

        fd_obj = ForestDisturbance(configuration_file_FD, geographic_boundaries, gdal_path,)

        fd_obj.run(FD_input_file_xml)
        
    logging.info("All outputs have been saved into: " + output_folder + "\n")
    logging.info("BIOMASS L2 Processor ended: see the above log messages for more info.")

    # copy the log and the configuration file to output
    if main_input_struct.proc_flags.AGB:
        log_file_name_new = os.path.join(os.path.dirname(log_file_name), "AGB", "biomassL2.log")
        copyfile(log_file_name, log_file_name_new)

        conf_file_name_new = os.path.join(os.path.dirname(log_file_name), "AGB", "ConfigurationFile.xml")
        copyfile(configuration_file_AGB, conf_file_name_new)

    if main_input_struct.proc_flags.FH:
        log_file_name_new = os.path.join(os.path.dirname(log_file_name), "FH", "biomassL2.log")
        copyfile(log_file_name, log_file_name_new)

        conf_file_name_new = os.path.join(os.path.dirname(log_file_name), "FH", "ConfigurationFile.xml")
        copyfile(configuration_file_FH, conf_file_name_new)

    if main_input_struct.proc_flags.TOMO_FH:
        log_file_name_new = os.path.join(os.path.dirname(log_file_name), "TOMO_FH", "biomassL2.log")
        copyfile(log_file_name, log_file_name_new)

        conf_file_name_new = os.path.join(os.path.dirname(log_file_name), "TOMO_FH", "ConfigurationFile.xml")
        copyfile(configuration_file_TOMO_FH, conf_file_name_new)

    return True


def start_logging(output_folder, proc_flags, log_level):
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

    log_file_name = os.path.join(output_folder, "biomassL2.log")

    logging.basicConfig(
        handlers=[logging.FileHandler(log_file_name, mode="w", encoding="utf-8"), logging.StreamHandler(),],
        level=level_to_set,
        format="%(asctime)s - %(levelname)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    logging.getLogger("matplotlib.font_manager").disabled = True

    logging.info(" --BIOMASS L2 Processor-- ")
    logging.info("Following chains will be executed: ")
    if proc_flags.AGB:
        logging.info("AGB")
    if proc_flags.FD:
        logging.info("FD")
    if proc_flags.FH:
        logging.info("FH (Interferometric phase)")
    if proc_flags.TOMO_FH:
        logging.info("FH (Tomographic phase)")
    if proc_flags.TOMO:
        logging.info("TOMO (Tomographic cube generation)")

    logging.info(" \n")

    return log_file_name


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

        output_folder_sub_chain = os.path.join(output_folder, "AGB")
        os.makedirs(output_folder_sub_chain)
        AGB_input_file_xml = os.path.join(output_folder_sub_chain, "InputFile.xml")

        input_params_obj = input_params(
            "AGB",
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
        AGB_input_file_xml = ""

    if main_input_struct.proc_flags.FD:

        output_folder_sub_chain = os.path.join(output_folder, "FD")
        os.makedirs(output_folder_sub_chain)
        FD_input_file_xml = os.path.join(output_folder_sub_chain, "InputFile.xml")

        input_params_obj = input_params(
            "FD",
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
        FD_input_file_xml = ""

    if main_input_struct.proc_flags.FH:

        output_folder_sub_chain = os.path.join(output_folder, "FH")
        os.makedirs(output_folder_sub_chain)
        FH_input_file_xml = os.path.join(output_folder_sub_chain, "InputFile.xml")

        input_params_obj = input_params(
            "FH",
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
        FH_input_file_xml = ""

    if main_input_struct.proc_flags.TOMO_FH:

        output_folder_sub_chain = os.path.join(output_folder, "TOMO_FH")
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
                    "The stacks with less than #3 acquisitions cannot be used to perform TOMO FH estimation \
                and they will be removed: #{} stacks among #{} totaL, have been removed".format(
                        len_stacks_tomo, len_stacks_all
                    )
                )

            TOMO_FH_input_file_xml = os.path.join(output_folder_sub_chain, "InputFile.xml")

            input_params_obj = input_params(
                "TOMO_FH",
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
            TOMO_FH_input_file_xml = ""
            logging.warning("cannot find any stack with more than #2 acquisitions: TOMO FH estimation will be disabled")
            main_input_struct.proc_flags.TOMO_FH = False
    else:
        TOMO_FH_input_file_xml = ""

    if main_input_struct.proc_flags.TOMO:

        output_folder_sub_chain = os.path.join(output_folder, "TOMO")
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
                    "The stacks with less than #3 acquisitions cannot be used to perform TOMO FH estimation \
                and they will be removed: #{} stacks among #{} totaL, have been removed".format(
                        len_stacks_tomo, len_stacks_all
                    )
                )

            TOMO_input_file_xml = os.path.join(output_folder_sub_chain, "InputFile.xml")

            input_params_obj = input_params(
                "TOMO",
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
            TOMO_input_file_xml = ""
            logging.warning("cannot find any stack with more than #2 acquisitions: TOMO estimation will be disabled")
            main_input_struct.proc_flags.TOMO = False

    else:
        TOMO_input_file_xml = ""

    return (
        AGB_input_file_xml,
        FD_input_file_xml,
        FH_input_file_xml,
        TOMO_FH_input_file_xml,
        TOMO_input_file_xml,
    )


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

            point_dict = main_input_struct.geographic_boundary.point

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
        ) = decode_unique_acquisition_id_string(stack_composition[unique_stack_id][0], output_format="string")

        unique_merged_stack_id = (
            "GC_" + global_cycle_idx + "_RGSW_" + rg_swath_idx + "_RGSBSW_" + rg_sub_swath_idx + "_AZSW_" + az_swath_idx
        )

        if unique_merged_stack_id in stacks_to_merge_dict.keys():
            stacks_to_merge_dict[unique_merged_stack_id].append(unique_stack_id)
        else:
            stacks_to_merge_dict[unique_merged_stack_id] = [unique_stack_id]

    return stacks_to_merge_dict


def main():
    from docopt import docopt

    args = docopt(__doc__, version="0.0.1")

    # Input file:
    input_file_xml = args["INPUTFILEXML"]
    conf_folder = args["--conf"]

    # run processor:
    biomassL2_processor_run(input_file_xml, conf_folder)


if __name__ == "__main__":
    main()
