# python imports
import os
import logging

from shutil  import  copyfile
# biomassL2 processor imports
from biomassL2.utility_orchestrator import start_logging, data_select_by_date_and_boundaries, write_chains_input_file_main, \
                                           check_if_path_exists, format_folder_name, collect_stacks_to_be_merged
from biomassL2.IO_interfaces        import parse_biomassL2_main_input_file

# following imports are tested with "try" because different modalities of 
# software releases may lack of some libraries, for example when chains are 
# delivered separately
try:
    from biomassL2.main_AGB     import main_AGB
except:
    pass
try:
    from biomassL2.main_FH      import main_FH
except:
    pass
try:
    from biomassL2.main_TOMO_FH import main_TOMO_FH
except:
    pass
try:
    from biomassL2.main_FD      import main_FD    
except:
    pass
try:
    from biomassL2.main_TOMO    import main_TOMO_CUBE
except:
    pass         


# main orchestrator:
# 1) Sets enviriment and logging, parses the main input
# 2) Chooses data to be processed
# 3) Launches each activated chain
def biomassL2_processor_main( input_file_xml, gdal_path, gdal_enviroment_path, INSTALLATION_FOLDER):
    
    # Set the enviroment
    os.environ['GDAL_DATA'] = gdal_enviroment_path
    
    
    # read the main input file
    main_input_struct = parse_biomassL2_main_input_file( input_file_xml )
     
    # date string format is: BIOMASS_L2_YYYYMMDDTHHMMSS
    current_date_string = format_folder_name()

    output_folder = os.path.join( main_input_struct.output_folder, current_date_string )
    
    if os.path.exists( output_folder ):
        error_message = 'Output Folder '+output_folder+' already exists, exiting now.'
        raise RuntimeError ( error_message )
    os.makedirs( output_folder )   
    
    
    # start the orchestrator logging
    log_file_name = start_logging( output_folder , main_input_struct.proc_flags, 'DEBUG')
    
    logging.debug( 'Installation   folder is {}'.format( INSTALLATION_FOLDER             ))
    logging.info( ' Auxiliary data folder is {}'.format( main_input_struct.L1c_aux_data_repository ))
    logging.info ( 'Results will be saved into output folder {}'.format( output_folder ))
    
    
    # get all the data that matches input datations and ground boundaries
    # a dictionary with stack_id and scene oroduct folder names is retrived
    logging.info('Searching data from data set : '+main_input_struct.L1c_repository )
    logging.info('Research will be done according to user dates and boundaries....\n')
    try:
        
        stack_composition, geographic_boundaries,  geographic_boundaries_per_stack = data_select_by_date_and_boundaries( main_input_struct )
        if not len(stack_composition):
            logging.error( 'Cannot find any data which meets the conditions to run the processing: see the log warnings. \n' )
            logging.info( 'BIOMASS L2 Processor ended without any processing performed.' )
            return 0
            
    except Exception as e:
        logging.error( e, exc_info=True )
        raise 
    
    logging.info('Geographic Boundaries extracted [deg]:')
    logging.info('Latitude  min: '+str( geographic_boundaries.lat_min ) )
    logging.info('Latitude  max: '+str( geographic_boundaries.lat_max ) )
    logging.info('Longitude min: '+str( geographic_boundaries.lon_min ) )
    logging.info('Longitude max: '+str( geographic_boundaries.lon_max ) )
    
    # write the input files for the activated chains:
    AGB_input_file_xml, FD_input_file_xml, FH_input_file_xml, TOMO_FH_input_file_xml, TOMO_input_file_xml = \
    write_chains_input_file_main( output_folder, main_input_struct, stack_composition )
  
    
    # Execute all the activated chains, separately
    default_configuration_folder = os.path.join(INSTALLATION_FOLDER, 'conf')
    check_if_path_exists(default_configuration_folder, 'FOLDER' )
    
    if main_input_struct.proc_flags.FH or main_input_struct.proc_flags.TOMO_FH:
        stacks_to_merge_dict = collect_stacks_to_be_merged( stack_composition )

    # AGB
    try:
    
        if main_input_struct.proc_flags.AGB:
            configuration_file_xml = os.path.join( default_configuration_folder, 'ConfigurationFile_AGB.xml')
            main_AGB( AGB_input_file_xml, configuration_file_xml, geographic_boundaries, geographic_boundaries_per_stack, gdal_path )
    except Exception as e:
        logging.error( 'Orchestrator: error inside AGB chain: '+str(e), exc_info=True )
        raise
   
    
    # FD    
    try:
        if main_input_struct.proc_flags.FD:
            configuration_file_xml = os.path.join( default_configuration_folder, 'ConfigurationFile_FD.xml')
            main_FD( FD_input_file_xml, configuration_file_xml, geographic_boundaries, gdal_path )
            
    except Exception as e:
        logging.error( 'Orchestrator: error inside FD chain: '+str(e), exc_info=True )
        raise
    
    
    # FH     
    try:    
        if main_input_struct.proc_flags.FH:
        
            configuration_file_xml = os.path.join( default_configuration_folder, 'ConfigurationFile_FH.xml')
            main_FH( FH_input_file_xml, configuration_file_xml, stacks_to_merge_dict, geographic_boundaries, gdal_path )
    
    except Exception as e:
        logging.error( 'Orchestrator: error inside FH chain: '+str(e), exc_info=True )
        raise
    
    
    # TOMO FH
    try:                 
        if main_input_struct.proc_flags.TOMO_FH:
            
            configuration_file_xml = os.path.join( default_configuration_folder, 'ConfigurationFile_TOMO_FH.xml')
            main_TOMO_FH( TOMO_FH_input_file_xml, configuration_file_xml, stacks_to_merge_dict, geographic_boundaries, gdal_path )
            
    except Exception as e:
        logging.error( 'Orchestrator: error inside TOMO FH chain: '+str(e), exc_info=True )
        raise
    
    
    # TOMO     
    try:        
        if main_input_struct.proc_flags.TOMO:
            configuration_file_xml = os.path.join( default_configuration_folder, 'ConfigurationFile_TOMO.xml')
            
            main_TOMO_CUBE(TOMO_input_file_xml, configuration_file_xml, gdal_path)
                
    except Exception as e:
        logging.error( 'Orchestrator: error inside TOMO chain: '+str(e), exc_info=True )
        raise


    logging.info( 'All outputs have been saved into: '+output_folder+'\n' )
    logging.info( 'BIOMASS L2 Processor ended: see the above log messages for more info.' )

    # copy the log and the configuration file to output
    if main_input_struct.proc_flags.AGB:
        log_file_name_new = os.path.join( os.path.dirname( log_file_name ), 'AGB', 'biomassL2.log' )
        copyfile(log_file_name, log_file_name_new)
        
        conf_file_name_curr = os.path.join( default_configuration_folder,            'ConfigurationFile_AGB.xml' )
        conf_file_name_new  = os.path.join( os.path.dirname( log_file_name ), 'AGB', 'ConfigurationFile.xml' )
        copyfile(conf_file_name_curr, conf_file_name_new)
        
    if main_input_struct.proc_flags.FD:  
        log_file_name_new = os.path.join( os.path.dirname( log_file_name ), 'FD',  'biomassL2.log' )
        copyfile(log_file_name, log_file_name_new)
        
        conf_file_name_curr = os.path.join( default_configuration_folder,            'ConfigurationFile_FD.xml' )
        conf_file_name_new  = os.path.join( os.path.dirname( log_file_name ), 'FD', 'ConfigurationFile.xml' )
        copyfile(conf_file_name_curr, conf_file_name_new)
        
    if main_input_struct.proc_flags.FH:
        log_file_name_new = os.path.join( os.path.dirname( log_file_name ), 'FH',  'biomassL2.log' )
        copyfile(log_file_name, log_file_name_new)
        
        conf_file_name_curr = os.path.join( default_configuration_folder,            'ConfigurationFile_FH.xml' )
        conf_file_name_new  = os.path.join( os.path.dirname( log_file_name ), 'FH', 'ConfigurationFile.xml' )
        copyfile(conf_file_name_curr, conf_file_name_new)
        
    if main_input_struct.proc_flags.TOMO_FH:
        log_file_name_new = os.path.join( os.path.dirname( log_file_name ), 'TOMO_FH', 'biomassL2.log' )
        copyfile(log_file_name, log_file_name_new)  
        
        conf_file_name_curr = os.path.join( default_configuration_folder,            'ConfigurationFile_TOMO_FH.xml' )
        conf_file_name_new  = os.path.join( os.path.dirname( log_file_name ), 'TOMO_FH', 'ConfigurationFile.xml' )
        copyfile(conf_file_name_curr, conf_file_name_new)
        
    if main_input_struct.proc_flags.TOMO:
        log_file_name_new = os.path.join( os.path.dirname( log_file_name ), 'TOMO', 'biomassL2.log' )
        copyfile(log_file_name, log_file_name_new)
        
        conf_file_name_curr = os.path.join( default_configuration_folder,            'ConfigurationFile_TOMO.xml' )
        conf_file_name_new  = os.path.join( os.path.dirname( log_file_name ), 'TOMO', 'ConfigurationFile.xml' )
        copyfile(conf_file_name_curr, conf_file_name_new)
        
        
    return True