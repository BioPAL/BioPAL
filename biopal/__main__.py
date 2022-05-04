"""
biopal processor

Usage:
  biopal [--conf=CONFFOLDER] INPUTFILEXML
  
Arguments:
  INPUTFILEXML       Path of the xml input file
  
Options:
  --conf=CONFFOLDER  Set the directory containing the xml configuration file
  -h --help          Show this screen
  --version          Show version
  
Details:
    processor can execute a single chain for each run;
    possible chains are AGB, FH, FD, TOMO_FH;
    the chain is set from INPUTFILEXML -> L2_product section

    If CONFFOLDER is not specified, default configurations are used
    
    Input and configuration can be generated with "biopal-quickstart"
"""

import os
import pkg_resources
import logging
import importlib
from biopal import __version__
from shutil import copyfile
from biopal.utility.utility_functions import (
    start_logging,
    check_if_path_exists,
    set_gdal_paths,
    format_folder_name,
)
from biopal.io.xml_io import parse_input_file, parse_configuration_file
from biopal.dataset_query.dataset_query import dataset_query

if not importlib.util.find_spec("biopal.agb.main_AGB") is None:
    from biopal.agb.main_AGB import AboveGroundBiomass
if not importlib.util.find_spec("biopal.fh.main_FH") is None:
    from biopal.fh.main_FH import ForestHeight
if not importlib.util.find_spec("biopal.tomo_fh.main_TOMO_FH") is None:
    from biopal.tomo_fh.main_TOMO_FH import TomoForestHeight
if not importlib.util.find_spec("biopal.fd.main_FD") is None:
    from biopal.fd.main_FD import ForestDisturbance
if not importlib.util.find_spec("biopal.fnf.main_FNF") is None:
    from biopal.fnf.main_FNF import ForestNonForestMap


class InvalidInputError(ValueError):
    pass


# main biopal
def biomassL2_processor_run(input_file_processor_xml, conf_folder=None):

    if conf_folder is None:
        configuration_file = pkg_resources.resource_filename("biopal", "_package_data/conf/Configuration_File.xml")
    else:
        default_configuration_folder = conf_folder
        check_if_path_exists(default_configuration_folder, "FOLDER")
        configuration_file = os.path.join(default_configuration_folder, "Configuration_File.xml")

    # parse the input file
    input_params_obj = parse_input_file(input_file_processor_xml)
    if (
        input_params_obj.L2_product is None
        or input_params_obj.output_specification is None
        or input_params_obj.dataset_query is None
    ):
        error_message = [
            "Main BioPAL APP input file should contain at least"
            ' "L2_product", "output_specification" and "dataset_query" sections'
        ]
        raise InvalidInputError(error_message)

    # date string format is: BIOMASS_L2_YYYYMMDDTHHMMSS
    current_date_string = format_folder_name()

    # add the date string sub-folder to output_folder
    output_folder = os.path.join(input_params_obj.output_specification.output_folder, current_date_string)
    # update the input params object output folder
    input_params_obj.output_specification.output_folder = output_folder

    if os.path.exists(output_folder):
        error_message = "Output Folder " + output_folder + " already exists, exiting now."
        raise RuntimeError(error_message)
    os.makedirs(output_folder)

    # start the main logging
    log_file_name = start_logging(output_folder, input_params_obj.L2_product, "DEBUG")

    logging.debug("Configuration file is {}".format(configuration_file))
    logging.info("Auxiliary data folder is {}".format(input_params_obj.dataset_query.L1C_aux_repository))
    logging.info("Results will be saved into output folder {}".format(output_folder))

    # get all the data that matches input datations and ground boundaries
    # a dictionary with stack_id and scene product folder names is retrived
    logging.info("Searching data from data set : " + input_params_obj.dataset_query.L1C_repository)
    logging.info("Research will be done according to user dates and boundaries....\n")
    try:

        dataset_query_obj = dataset_query()
        input_file_for_stack_based = dataset_query_obj.run(input_params_obj)

    except Exception as e:
        logging.error(e, exc_info=True)
        raise

    # set the GDAL environment
    conf_params_obj = parse_configuration_file(configuration_file)
    set_gdal_paths(conf_params_obj.gdal.gdal_path, conf_params_obj.gdal.gdal_environment_path)

    # Execute the activated chain
    if input_params_obj.L2_product == "AGB":
        chain_obj = AboveGroundBiomass(configuration_file)

    if input_params_obj.L2_product == "FH":
        chain_obj = ForestHeight(configuration_file)

    if input_params_obj.L2_product == "TOMO_FH":
        chain_obj = TomoForestHeight(configuration_file)

    if input_params_obj.L2_product == "FD":
        chain_obj = ForestDisturbance(configuration_file)
        
    if input_params_obj.L2_product == "FNF":
        chain_obj = ForestNonForestMap(configuration_file)

    chain_obj.run(input_file_for_stack_based)

    # copy the configuration file to output
    input_file_name_out = os.path.join(os.path.dirname(log_file_name), "Input_File.xml")
    copyfile(input_file_processor_xml, input_file_name_out)

    conf_file_name_out = os.path.join(os.path.dirname(log_file_name), "Configuration_File.xml")
    copyfile(configuration_file, conf_file_name_out)

    logging.info("Outputs have been saved into: " + output_folder + "\n")
    logging.info("BIOMASS L2 Processor ended: see the above log messages for more info.")

    return True


def main():
    from docopt import docopt

    args = docopt(__doc__, version=__version__)

    # Input file:
    input_file_processor_xml = args["INPUTFILEXML"]
    conf_folder = args["--conf"]

    # run processor:
    biomassL2_processor_run(input_file_processor_xml, conf_folder)


if __name__ == "__main__":
    main()
