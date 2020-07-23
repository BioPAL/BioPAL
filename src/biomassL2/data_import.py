import os
import logging
import numpy as np
from   scipy.interpolate           import interp1d
from   biomassL2.IO_interfaces     import raster_info
from   arepytools.io.productfolder import ProductFolder


def get_data_time_stamp( folder, pf_name ):
    # reads a data:
    # it is supposed to contain one or more polarizations (the "SwathInfo" is read to retrive it)
    # it returns a dictionary with keys = polarizations 
    # it returns also the dimensions of data
    
    data_pf_name = os.path.join( folder, pf_name)

    pf = ProductFolder(data_pf_name, 'r')

    # prepare the metadata elements
    data_channel_obj    = pf.get_channel( 0 )
    metadata_obj        = data_channel_obj.metadata
    metadatachannel_obj = metadata_obj.get_metadata_channels(0)
    
    # Raster Info
    ri = metadatachannel_obj.get_element('RasterInfo')

    lines_start_utc        = str(ri.lines_start)
    
    return  lines_start_utc



def read_data( folder, pf_name ):
    # reads a data:
    # it is supposed to contain one or more polarizations (the "SwathInfo" is read to retrive it)
    # it returns a dictionary with keys = polarizations 
    # it returns also the dimensions of data
    
    data_pf_name = os.path.join( folder, pf_name)

    pf = ProductFolder(data_pf_name, 'r')
    number_of_pols = pf.get_number_channels()
    data_read = {}
    polid_found = []
    for pol_channel_idx in range( number_of_pols ):
        
        # prepare the metadata elements
        data_channel_obj    = pf.get_channel( pol_channel_idx )
        metadata_obj        = data_channel_obj.metadata
        metadatachannel_obj = metadata_obj.get_metadata_channels(0)
        
        # get the ID of the master acquisition:
        di = metadatachannel_obj.get_element('DataSetInfo')
        if not di:
            raise ValueError('data product folder should contain the DataSetInfo to retrive the MASTER ID')
        if di.description.find( 'Master_swath_' ) != 0:
            raise ValueError('DataSetInfo description not recognized: it should be a string as "Master_swath_IdOfTheMaster"' )
        master_id = di.description[13:]
        
        # Raster Info
        ri = metadatachannel_obj.get_element('RasterInfo')
        num_samples            = ri.samples
        num_lines              = ri.lines
        pixel_spacing_slant_rg = ri.samples_step
        pixel_spacing_az       = ri.lines_step
        lines_start_utc        = str(ri.lines_start)
        
        # SwathInfo
        si = metadatachannel_obj.get_element('SwathInfo')
        if not si:
            raise ValueError('data product folder should contain the SwathInfo to retrive the polarization')
        pol_id = si.polarization.name
        
        polid_found.append(pol_id)
        # Sampling constants
        sc = metadatachannel_obj.get_element('SamplingConstants')
        resolution_m_slant_rg = sc.brg_hz
        resolution_m_az       = sc.baz_hz
        
        # hv and vh data are mean togheter, ew save only a vh polarization, that will be a "vh_used = (vh+hv)/2"
        if pol_id == 'hv' or pol_id == 'vh':
            if 'vh' in data_read.keys():
                # data (vh or hv) already saved in the dict, add the other data
                data_read[ 'vh' ] = ( data_read[ 'vh' ] + pf.read_data(pol_channel_idx).transpose() )/np.sqrt(2)
            else:
                # nor vh nor vv have been saved to dict yet, add first one
                data_read[ 'vh' ] = pf.read_data(pol_channel_idx).transpose() 
        else:
        
            data_read[ pol_id ] = pf.read_data(pol_channel_idx).transpose()  

    if len(polid_found)<4:
        raise ValueError('Input data stack {} should contain #4 polarizations, hh, hv, vh, vv, only {} found '.format( pf_name, len(polid_found) ) )
    elif not 'hh' in polid_found or not 'hv' in polid_found or not 'vh' in polid_found or not 'vv' in polid_found:
        raise ValueError('Input data stack {} should contain #4 polarizations, hh, hv, vh, vv, only {} found '.format( pf_name, len(polid_found) ) )
                                                                                                                                                          
    return data_read, num_samples, num_lines, pixel_spacing_slant_rg, pixel_spacing_az, resolution_m_slant_rg, resolution_m_az, master_id, lines_start_utc


def data_oversample( data, oversampling_factor, raster_info_obj ):
    
    rg_ratio = np.floor( raster_info_obj.resolution_m_slant_rg / raster_info_obj.pixel_spacing_slant_rg )
    az_ratio = np.floor( raster_info_obj.resolution_m_az / raster_info_obj.pixel_spacing_az )
        
    rg_oversampling_flag = False
    az_oversampling_flag = False
    
    if rg_ratio<2:
        rg_oversampling_flag       = True
        num_samples_out            = raster_info_obj.num_samples * oversampling_factor
    else:
        num_samples_out            = raster_info_obj.num_samples
        pixel_spacing_slant_rg_out = raster_info_obj.pixel_spacing_slant_rg
        
    if az_ratio<2:
        az_oversampling_flag = True
        num_lines_out        = raster_info_obj.num_lines * oversampling_factor
    else:
        num_lines_out        = raster_info_obj.num_lines
        pixel_spacing_az_out = raster_info_obj.pixel_spacing_az
    
    
    if rg_oversampling_flag or az_oversampling_flag:
        logging.info('    Oversampling needed:')

        if not type(data) is dict:
            # slope and reference_heights cases  (input data is value matrix)
            
            if rg_oversampling_flag:
                
                logging.info('        range   oversampling of auxiliary data...')
                
                data, pixel_spacing_slant_rg_out = data_oversample_core( data, 0, raster_info_obj.pixel_spacing_slant_rg, oversampling_factor )
                
            if az_oversampling_flag:
                
                logging.info('        azimuth   oversampling of auxiliary data...')

                data, pixel_spacing_az_out = data_oversample_core( data, 1, raster_info_obj.pixel_spacing_az, oversampling_factor )
        
        else:
            
            for data_key, data_extracted in data.items():
                
                if not type(data_extracted) is dict:
                    

                    # KZ, off_nadir and ECEFGRID case (input data is a dict of values)
                    if rg_oversampling_flag:
                        
                        logging.info('        range   oversampling of '+data_key+'...')
                        
                        data[data_key], pixel_spacing_slant_rg_out = data_oversample_core( data_extracted, 0, raster_info_obj.pixel_spacing_slant_rg, oversampling_factor )
                
                    if az_oversampling_flag:
                        logging.info('        azimuth oversampling of '+data_key+'...')
                        
                        data_extracted = data[data_key]
                        data[data_key], pixel_spacing_az_out = data_oversample_core( data_extracted, 1, raster_info_obj.pixel_spacing_az, oversampling_factor )

                else:
                    
                    # Beta0 case (input data is a dict of dict with values)
                    for pol_key, data_pol in data_extracted.items():
                        
                        if rg_oversampling_flag:
                            
                            logging.info('        range   oversampling of '+data_key+' , polarization '+pol_key+'...')

                            data[data_key][pol_key], pixel_spacing_slant_rg_out = data_oversample_core( data_pol, 0, raster_info_obj.pixel_spacing_slant_rg, oversampling_factor )
                
                        if az_oversampling_flag:
                            logging.info('        azimuth oversampling of '+data_key+' , polarization '+pol_key+'...')
                            
                            data_pol = data[data_key][pol_key] 
                            data[data_key][pol_key], pixel_spacing_az_out = data_oversample_core( data_pol, 1, raster_info_obj.pixel_spacing_az, oversampling_factor )

            
        logging.info('    oversampling done.\n')
                
    return data, num_samples_out, pixel_spacing_slant_rg_out, num_lines_out, pixel_spacing_az_out



def data_oversample_core( data, axis_index, pixel_spacing_in, oversampling_factor ):
    
    # original size of the axis to be interpolated
    axis_len_in = data.shape[ axis_index ]
    
    # pixel spacing after the interpolation
    pixel_spacing_out = pixel_spacing_in / oversampling_factor
    
    # input original axis
    max_val_in = axis_len_in*pixel_spacing_in
    ax_in      = np.arange( 0, max_val_in, pixel_spacing_in)
    
    # output interpolated axis 
    max_val_out = axis_len_in*oversampling_factor*pixel_spacing_out
    ax_out      = np.arange( 0, max_val_out, pixel_spacing_out)
    
    # interpolation of data along axis_index
    interp_fun = interp1d (ax_in,   data, axis=axis_index, bounds_error=0)
    data       = interp_fun( ax_out )
    
    return data, pixel_spacing_out


def read_auxiliary_single_channel( folder, pf_name ):
    # reads Incidence_Angle and Reference_height auxiliary data:
    # it returns a numpy matrix, no dictionaries in this case
    
    data_pf_name = os.path.join( folder, pf_name)
    
    if os.path.exists( data_pf_name ):
        pf = ProductFolder(data_pf_name, 'r')
        number_of_channels = pf.get_number_channels()
        if number_of_channels>1:
            raise ValueError('Input auxiliary data is supposed to have just one channel, and not # {}'.format(number_of_channels)) 
        
        aux_read = pf.read_data( 0 ).transpose()   
    else:
        aux_read = None
        logging.warning( 'Path '+data_pf_name+' does not exist.' )
        
    return aux_read


def read_auxiliary_multi_channels( folder, pf_name, valid_acq_id_to_read=None, read_raster_info=False ):
    # reads a KZ product:
    # it is supposed to be a pf containing "N" channels, with "N" the number of acquisitions in a stack
    # the acquisition name is read from the SwathInfo "Swath" field
    # it returns a dictionary with keys = acquisition_id ( which is the "Swath")
    
    data_pf_name = os.path.join( folder, pf_name)

    if os.path.exists( data_pf_name ):
        pf = ProductFolder(data_pf_name, 'r')
        number_of_acq = pf.get_number_channels()

        data_read = {}
        for channel_idx in range( number_of_acq ):
            
            # prepare the metadata elements
            data_channel_obj    = pf.get_channel( channel_idx )
            metadata_obj        = data_channel_obj.metadata
            metadatachannel_obj = metadata_obj.get_metadata_channels(0)
            
            # SwathInfo
            si = metadatachannel_obj.get_element('SwathInfo')
            if not si:
                raise ValueError('Input KZ and off_nadir should contain the SwathInfo to retrive the Swath ID')
             
            if valid_acq_id_to_read is None or ( si.swath in valid_acq_id_to_read ):

                data_read[ si.swath ] = pf.read_data(channel_idx).transpose()  
        
        # Raster Info
        ri = metadatachannel_obj.get_element('RasterInfo')
        num_samples            = ri.samples
        num_lines              = ri.lines
        pixel_spacing_slant_rg = ri.samples_step
        pixel_spacing_az       = ri.lines_step      
        
        raster_info_obj = raster_info ( num_samples,
                                num_lines,
                                pixel_spacing_slant_rg,
                                pixel_spacing_az,
                                None,
                                None,
                                None
                                )
    else:
        data_read       = None
        raster_info_obj = None
        logging.warning( 'Path '+data_pf_name+' does not exist.' )
      
    if read_raster_info:
        return data_read, raster_info_obj
    else:
        return data_read


def read_ecef_grid( folder, pf_name ):
    # reads an ECEFGRID:
    # it is supposed to be a pf containing exactly 3 channels, with X, Y and Z coordinates
    # the X,Y or Z is got from the DataSetInfo Description field which is supposed
    # to be exactly a string like :Auxiliary data: X ECEF GRID [m] (example for X)
    # it returns a dictionary with keys = coordinate_id (X, Y or Z) 
    
    data_pf_name = os.path.join( folder, pf_name)
    
    if os.path.exists( data_pf_name ):
        pf = ProductFolder(data_pf_name, 'r')
        number_of_coords = pf.get_number_channels()
        if not number_of_coords==3:
            raise ValueError('Input ECEF GRID should contain #3 channels with X,Y and Z coordinates: #{}  channels have been found.'.format(number_of_coords))
                             
        coordinates_read = {}
        for coord_channel_idx in range( number_of_coords ):
            
            # prepare the metadata elements
            data_channel_obj    = pf.get_channel( coord_channel_idx )
            metadata_obj        = data_channel_obj.metadata
            metadatachannel_obj = metadata_obj.get_metadata_channels(0)
            
            # DataSetInfo
            di = metadatachannel_obj.get_element('DataSetInfo')
            if not di:
                raise ValueError('Input ECEF GRID should contain the DataSetInfo to retrive the Description')
            coord_id = di.description[16]
            if not coord_id=='X' and not coord_id=='Y' and not coord_id=='Z':
                raise ValueError('Cannot retrive coordinate name from DataSetInfo description: description should be a string as: "Auxiliary data: X ECEF GRID [m]", instead it is: "'+di.description+'"')
            
            coordinates_read[ coord_id ] = pf.read_data(coord_channel_idx).transpose()  
            
    else:
        coordinates_read = None
        logging.warning( 'Path '+data_pf_name+' does not exist.' ) 
        
    return coordinates_read