# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 14:48:46 2020

@author: macie
"""

import numpy as np
import scipy as sp
import logging
import os

from biopal.agb.processing_AGB import (
    mean_on_rois,
    check_intersection,
    interp2d_wrapper,
    merge_agb_intermediate
)


# %% transform functions
def transform_function(in_data,interval,kind,do_forward=True):
    out_data = np.copy(in_data)
    out_data[(out_data<interval[0]) | (out_data>interval[1])] = np.nan
    # note: no check of data and kind is done here, 
    # it is assumed that the inputs are correct
    if kind=='db':
        if do_forward:
            out_data = 10*np.log10(out_data)
        else:
            out_data = 10**(0.1*in_data)
    elif kind=='-db':
        if do_forward:
            out_data = -10*np.log10(in_data)
        else:
            out_data = 10**(-0.1*in_data)
    elif kind=='-2db':
        if do_forward:
            out_data = -10*np.log10(2*in_data)
        else:
            out_data = 10**(-0.1*in_data)/2
    elif kind=='cosdb':
        if do_forward:
            out_data = 10*np.log10(np.cos(in_data))
        else:
            out_data = np.arccos(10**(0.1*in_data))
    else:
        out_data = np.copy(in_data)
    return out_data


# %% fitting to data table
def fit_formula_to_table_data(original_formula,
                              observable_table,
                              observable_table_column_names,
                              parameter_position_table,
                              parameter_position_table_column_names,
                              individual_parameter_min_max_tables):
    
    # it is required that observable_table contains no nans in columns that double with parameter_position_table
    # (in terms of columns names)
    # in columns that double, those rows that have observable values are treated as calibration data
    # no other unnecessary columns are allowed in observable_table and parameter_position_table (i.e., with names
    # that do not occur in original_formula)
    
    # convert the columnwise indices in "parameter_position_table" to unique indices
    unique_index_table,columnwise_to_unique_index_lut = \
        regularise_indices(parameter_position_table)
        
    # number of unique parameters to estimate
    length_of_p_vector = columnwise_to_unique_index_lut.shape[0]
    
    # create a table of min and max parameter values (requires looping through parameter tables and extracting relevant rows)
    # this allows a possible flexible setting of intervals in the future
    p_min_max_table = np.nan * np.zeros((length_of_p_vector,2))
    for parameter_idx,individual_parameter_min_max_table in enumerate(individual_parameter_min_max_tables):
        current_rows_in_lut = columnwise_to_unique_index_lut[:,1]==parameter_idx
        current_positions_in_individual_parameter_table = np.int32(columnwise_to_unique_index_lut[current_rows_in_lut,2])
        p_min_max_table[current_rows_in_lut,:] = \
            individual_parameter_min_max_table[current_positions_in_individual_parameter_table,:]
    # extract min, max, intiial values
    p_lower = p_min_max_table[:,0]
    p_upper = p_min_max_table[:,1]
    
    # the name of the table that is used in the original_formula string 
    # (could be anything as long as the same string is used in swap_names_in_original_formula and cost_function)
    table_name_in_converted_formula = "current_data_table"
    
    
    # find rows for which all observable data exist
    rows_with_all_observables = np.all(~np.isnan(observable_table),axis=1)
    # in this case, the formula should use the first occurence of the same quantity, if it is observed in both "observables" and "parameters"
    converted_calibration_formula = swap_names_and_merge_formula(original_formula,
                                                observable_table_column_names,
                                                parameter_position_table_column_names,
                                                table_name_in_converted_formula,
                                                use_observable_if_repeated_and_available = True)
    
    if np.all(rows_with_all_observables):
        # if all rows have all data, the total cost function uses only the calibration formula
        cost_function_arguments = ([converted_calibration_formula],
                [observable_table],
                [unique_index_table],
                table_name_in_converted_formula,
                lambda x:parameter_transfer_function(x,p_lower,p_upper,False))
    else:
        # add an estimation ssd
        converted_estimation_formula = swap_names_and_merge_formula(original_formula,
                                            observable_table_column_names,
                                            parameter_position_table_column_names,
                                            table_name_in_converted_formula,
                                            use_observable_if_repeated_and_available = False)
        cost_function_arguments = ([converted_calibration_formula,converted_estimation_formula],
                [observable_table[rows_with_all_observables,:],observable_table],
                [unique_index_table[rows_with_all_observables,:],unique_index_table],
                table_name_in_converted_formula,
                lambda x:parameter_transfer_function(x,p_lower,p_upper,False))
    
        
    
    # iterate a few times with different initial values in case some initial value set fails
    for counter in range(10):
        
        # creating initial values by randomising
        p_initial = p_lower+np.random.rand(length_of_p_vector)*(p_upper-p_lower)
        # converting to x
        x_initial = parameter_transfer_function(p_initial,p_lower,p_upper,True)
        
        # fit the model
        fitted_model = sp.optimize.minimize(cost_function,x_initial,cost_function_arguments,method='BFGS')
        print(fitted_model.message)
        if fitted_model.success:
            p_estimated = parameter_transfer_function(fitted_model.x,p_lower,p_upper,False)
            break
        else:
            p_estimated = np.nan * np.zeros(length_of_p_vector)
    return (
        np.column_stack((columnwise_to_unique_index_lut,p_initial,p_estimated)),
        p_estimated[unique_index_table],
        cost_function(parameter_transfer_function(p_estimated,p_lower,p_upper,True),*cost_function_arguments))

def match_string_lists(ref_string_list,test_string_list):
    is_in = np.zeros((len(ref_string_list),len(test_string_list)))
    for ref_idx,current_ref_string in enumerate(ref_string_list):
        for test_idx,current_test_string in enumerate(test_string_list):
            is_in[ref_idx,test_idx] = current_ref_string.find(current_test_string)
    return is_in

# %% function for calculating a given statistic for polygons of arbitrary shape  
def stats_on_polygons(data,pixel_axis_east,pixel_axis_north,reference_polygons,method):
    # here, create a function that calculates statistic in method on CAL data polygons
    #
    print('not implemented')
    return []
# %% function for calculating a given statistic on all samples on a grid and all polygons
def stats_on_all_samples(data,pixel_axis_east,pixel_axis_north,sample_axis_east,sample_size_east,sample_axis_north,sample_size_north,polygons,method):
    stats = mean_on_rois(data,pixel_axis_east,pixel_axis_north,sample_axis_east,sample_size_east,sample_axis_north,sample_size_north,method)
    if polygons:
        stats_polygons = stats_on_polygons(data,pixel_axis_east,pixel_axis_north,polygons,method)
        stats = np.concatenate((stats,stats_polygons))
    return stats

# %%
# function for creating list of format strings and a header for subsequent saving of tables into text files
def get_fmt_and_header(column_names,all_column_groups,all_data_types,delimiter='\t',precision=2,column_width=10):
    out_format = []
    out_header = []
    for curr_column in column_names:
        for curr_column_group,curr_data_type in zip(all_column_groups,all_data_types):
            if curr_column in curr_column_group:
                if curr_data_type == 'd':
                    out_format.append('%s%dd' % (r'%',column_width))
                elif curr_data_type == 'f':
                    out_format.append('%s%d.%df' % (r'%',column_width,precision))
                break
        out_header.append('%s%ds' % (r'%',column_width) % curr_column)
    return out_format,delimiter.join(out_header)

#%% swap variable names in formulas to slices of an array
def swap_names_and_merge_formula(original_formulas,observable_names,parameter_names,new_table_name,use_observable_if_repeated_and_available = True):
    original_variable_names = observable_names + parameter_names
    unique_variable_names,name_counts = np.unique(np.array(original_variable_names),return_counts = True)
    new_formula = '0'
    for current_formula in original_formulas:
        for unique_variable_name,name_count in zip(unique_variable_names,name_counts):
            if name_count==1:
                position_in_variable_names_vector = np.where(np.array(original_variable_names)==unique_variable_name)[0][0]
            elif name_count==2:
                position_in_variable_names_vector = np.where(np.array(original_variable_names)==unique_variable_name)[0][np.int32(~use_observable_if_repeated_and_available)]
            current_formula = current_formula.replace(unique_variable_name,new_table_name + ('[:,%d]' % (position_in_variable_names_vector)))
        new_formula = new_formula + ' + (%s)**2' % (current_formula)
    return new_formula.replace('0 + ','')
    
# %% cost function from formulas and tables
def cost_function(x_vector,converted_formulas,observable_tables,index_tables,name_of_table_in_converted_formula,transfer_function):
    # number of independent formulas
    number_of_formulas = len(converted_formulas)
    # convert unconstrained x vector to constrained p vector
    p_vector = transfer_function(x_vector)
    # allocate total cost
    total_cost = 0
    # loop through different formulas and associated tables
    for observable_table,index_table,converted_formula in zip(observable_tables,index_tables,converted_formulas):
        table_in_converted_formula = np.column_stack((observable_table,p_vector[index_table]))
        # evaluate the cost function formula, average across different measurements and add to the total cost
        total_cost += np.mean(eval(converted_formula,{name_of_table_in_converted_formula:table_in_converted_formula}))
        # print(eval(converted_formula,{name_of_table_in_converted_formula:table_in_converted_formula}))
    # return weighed with the number of formulas    
    return np.sqrt(total_cost/number_of_formulas)   

#%% function for converting columnwise indices (which are repeated within the same column if they represent identical values,
# but which may be repeated across different columns without meaning that they represent identical values)
# to unique indices (which are only repeated within the same table if they are to have identical values)
def regularise_indices(columnwise_index_table):
    
    offset = 0
    unique_index_table = [] # position in a single x-vector
    columnwise_to_unique_index_lut = [] # lut for converting between parameter id and column and parameter position in x
    for column_idx,parameter_column in enumerate(columnwise_index_table.transpose()):
        # add -1 at the beginning to avoid vectors with one element (which will not work with interp1d)
        old_indices = np.concatenate((-1*np.ones(1),np.unique(parameter_column)))
        # new indices is a simple sequence from 0 to number of parameters-1 + offset due to previous parameters
        new_indices = np.arange(len(old_indices))-1+offset
        # convert parameter indices and add to the list
        unique_index_table.append(sp.interpolate.interp1d(old_indices,new_indices,kind='nearest')(parameter_column))
        # save the lut, removing the first, unnecessary element
        columnwise_to_unique_index_lut.append(
            np.column_stack((new_indices[old_indices>-1],
                             column_idx+np.zeros(len(old_indices[old_indices>-1])),
                             old_indices[old_indices>-1],
                             )))
        # update offset based on current parameter column
        offset = np.max(new_indices)+1
    # convert the list of vectors to an array
    unique_index_table = np.int32(np.column_stack(unique_index_table))
    # stack all luts to one lut
    columnwise_to_unique_index_lut = np.row_stack(columnwise_to_unique_index_lut)
    # # length of beta vector
    # length_of_p_vector = columnwise_to_unique_index_lut.shape[0]
    return unique_index_table,columnwise_to_unique_index_lut

# %% define parameter transfer functions
def parameter_transfer_function(in_vector,p_lower,p_upper,in_vector_is_p=False):
    # note: no check of x_vector, p_upper, p_lower is done here, 
    # it is assumed that the inputs are correct
    if not in_vector_is_p:
        return p_lower+(p_upper-p_lower)*np.sin(in_vector)**2
    else:
        return np.arcsin(np.sqrt((in_vector-p_lower)/(p_upper-p_lower)))


# %% reading tables
def read_tables(
        block_extents, # extent of the current area for which the table is created
        pixel_axis_east, # east north axes onto which data are interpolated
        pixel_axis_north,
        sample_axis_east, # east north axes of the regular sampling grid
        sample_axis_north,
        sample_size_east, # east north extents of the samples
        sample_size_north,
        additional_sampling_polygons, # additional arbitrarily shaped polygons
        stack_in_block, # flags whether each stack is in current block
        stack_info_table, # info table with stack properties (stack id, headings, etc., defining the acquisition parameters)
        stack_info_table_column_names, # column names for the abovementioned table
        forest_class_boundaries, 
        forest_class_paths, 
        observable_names, # observable names in formula
        observable_is_stacked, # flags whether the observable comes from SAR and should be interpreted with stack_info_table or should be replicated across stacks
        observable_is_required, # flags whether the observable must exist in each table row
        observable_source_paths, # paths to equi7 tiles or files to read for each observable
        observable_source_bands, # which band in the tiff file to read (1 = first band)
        observable_transforms, # transform function to apply to observable
        observable_averaging_methods, # averaging method (most commonly 'mean', but could be other if required (e.g., for slope aspect angle))
        observable_ranges, # permitted ranges, outside those the observable is set to nan
        parameter_names, # parameter names in formula
        parameter_limits, # permissible parameter intervals
        parameter_variabilities, # parameter variabilities across all dimensions
        number_of_subsets, # number of subsets to use (used to allocate columns in parameter tables)
        ):

    
    # derived parameters
    number_of_stacks = stack_info_table.shape[0]
    stack_id_vector = stack_info_table[:,0]
    number_of_samples = len(sample_axis_east)*len(sample_axis_north) + len(additional_sampling_polygons)
    number_of_parameters = len(parameter_names)
    number_of_observables = len(observable_names)
    pixel_mesh_east, pixel_mesh_north = np.meshgrid( pixel_axis_east, 
                                                            pixel_axis_north)
            
    
    # defining names for identifiers (sampleID & forest class ID and then all columns from stack info table)
    identifier_names = ['sample_id','forest_class_id'] + stack_info_table_column_names
    number_of_identifiers = len(identifier_names)
    identifier_table = np.nan * np.zeros (( 
        number_of_samples * number_of_stacks,
        number_of_identifiers
        ))
    
    
        
    # filling out the first column with sample IDs
    identifier_table[:,0] = np.kron(np.arange(number_of_samples),np.ones(number_of_stacks))
    
    # filling out columns 3-8 with stack IDs and corresponding other identifications from stack_info_table
    identifier_table[:,2] = np.kron(np.ones(number_of_samples),stack_id_vector)
    for id_idx in range(5):
        identifier_table[:,3+id_idx] = sp.interpolate.interp1d(stack_id_vector,stack_info_table[:,1+id_idx],kind='nearest')(identifier_table[:,2])
    
    
    ### READING AND SAMPLING FOREST CLASS DATA
    logging.info('reading & sampling forest class data')
    # get equi7 information
    # equi7_grid_name = os.listdir(os.path.join(stack_paths[0], 'sigma0_hh'))[0]
    # loading all forest class maps that fall inside a block: note, more than one equi7 tile can fall inside
    # loaded masks needs to be re-interpolated over the pixels grid (pixel_axis_east, pixel_axis_north)
    forest_class_map_interp = np.zeros(
        (len(pixel_axis_north), len(pixel_axis_east)), dtype='float'
    )

    counter_forest_class_maps = 0
    for forest_class_map_idx, forest_class_map_path in enumerate(forest_class_paths):

        forest_class_map_boundaries = forest_class_boundaries[forest_class_map_idx]

        forest_class_map_is_inside = check_intersection(
            forest_class_map_boundaries[0],
            forest_class_map_boundaries[1],
            forest_class_map_boundaries[2],
            forest_class_map_boundaries[3],
            block_extents[0],
            block_extents[1],
            block_extents[3],
            block_extents[2],
        )
        if forest_class_map_is_inside:

            forest_class_map_interp_curr = np.round(
                interp2d_wrapper(
                    forest_class_map_path, 1, pixel_axis_east, pixel_axis_north, fill_value=float(0)
                )
            )

            # mean all the fnf tiles
            forest_class_map_interp = np.ceil(
                merge_agb_intermediate(
                    forest_class_map_interp, forest_class_map_interp_curr, method='nan_mean'
                )
            )
            
            # set all unrealistic values to 0 = non-forest
            forest_class_map_interp[(forest_class_map_interp <= 0) | np.isnan(forest_class_map_interp)] = 0
            
            counter_forest_class_maps = counter_forest_class_maps + 1

    if counter_forest_class_maps == 0:

        err_str = 'Cannot find any forest class map falling in current block coordinates.'
        logging.error(err_str)
        raise ImportError(err_str)
    
    # sampling forest class map
    temp_forest_class_vector = np.int32(np.round(stats_on_all_samples(
        forest_class_map_interp,
        pixel_mesh_east,
        pixel_mesh_north,
        sample_axis_east,
        sample_size_east,
        sample_axis_north,
        sample_size_north,
        additional_sampling_polygons,
        'mode',
    )))
    # repeating it across stacks and inserting into observable data table (note: forest class assumed constant across stacks)
    identifier_table[:,1] = np.kron(temp_forest_class_vector,np.ones(number_of_stacks))
    
    
    
    
    
    
    ### READING OBSERVABLE DATA
    # allocate observable table
    observable_table = np.nan * np.zeros((
        number_of_samples * number_of_stacks,
        number_of_observables))
    # cycle through observable sources in the stack
    for observable_idx in range(number_of_observables):
        
        # check if observable is stacked (if yes, cycle through stacks)
        if observable_is_stacked[observable_idx]:
                
            # cycle over stacks
            for stack_idx in range(number_of_stacks):
                
        
                # go ahead only if current stack is (at least partially) contained in the current parameter block:
                if stack_in_block[stack_idx]:
                    
                    # extracting various parts of the path
                    equi7_tiles_names = os.listdir(observable_source_paths[observable_idx][stack_idx])
                    current_head,equi7_grid_name = os.path.split(observable_source_paths[observable_idx][stack_idx])
                    current_head,source_name = os.path.split(current_head)
                    base_name = os.path.split(current_head)[1]
                    
                    
                    # cycle over all the equi7 tiles, interpolate over pixel grid and average
                    source_data_interp = np.NaN * np.zeros(
                        (len(pixel_axis_north), len(pixel_axis_east)), dtype='float'
                    )
                    for equi7_tile_name in equi7_tiles_names:
                        
                        
                        source_tiff_name = os.path.join(
                            observable_source_paths[observable_idx][stack_idx],
                            equi7_tile_name,
                            source_name 
                            + '_'
                            + base_name
                            + '_'
                            + equi7_grid_name[6:]
                            + '_'
                            + equi7_tile_name
                            + '.tif',
                        )

                        source_data_interp_curr = interp2d_wrapper(
                            source_tiff_name,
                            observable_source_bands[observable_idx],
                            pixel_axis_east,
                            pixel_axis_north,
                            fill_value=np.NaN,
                        )

                        source_data_interp = merge_agb_intermediate(
                            source_data_interp, source_data_interp_curr, method='nan_mean'
                        )

                    # masking the stack:
                    source_data_interp[forest_class_map_interp == 0] = np.NaN

                    logging.info("sampling and transforming data in file (band {}): \n\t'{}'\nto observable '{}' using transform '{}' and averaging method '{}' (mean value: {})".format(
                        observable_source_bands[observable_idx],
                        source_tiff_name,
                        observable_names[observable_idx],
                        observable_transforms[observable_idx],
                        observable_averaging_methods[observable_idx],
                        np.nanmean(source_data_interp)))
                    
                    # calculate sample statistics
                    temp_transformed_sampled_data = transform_function(
                        stats_on_all_samples(
                            source_data_interp,
                            pixel_mesh_east,
                            pixel_mesh_north,
                            sample_axis_east,
                            sample_size_east,
                            sample_axis_north,
                            sample_size_north,
                            additional_sampling_polygons, 
                            observable_averaging_methods[observable_idx],
                            ),
                        observable_ranges[observable_idx],
                        observable_transforms[observable_idx])
                    # find rows where the stack id matches the current stack id
                    current_rows = (identifier_table[:,2]==stack_id_vector[stack_idx])
                    # fill out the table
                    observable_table[current_rows,
                                          observable_idx] = temp_transformed_sampled_data

        
        # otherwise, do not cycle through stacks  
        # simply read the specified files or list of files
        else:
    
            source_data_interp = np.nan * np.zeros(
                (len(pixel_axis_north), len(pixel_axis_east)), dtype='float'
            )
            
        
            for file_idx,file_path in enumerate(observable_source_paths[observable_idx]):

                source_data_interp_curr = interp2d_wrapper(
                    file_path, observable_source_bands[observable_idx], pixel_axis_east, pixel_axis_north, fill_value=np.NaN
                )
                # mean all the equi7 tiles
                source_data_interp = merge_agb_intermediate(
                    source_data_interp, source_data_interp_curr, method='nan_mean'
                )
                
                logging.info("sampling and transforming data in file (band {})): \n\t'{}'\nto observable '{}' using transform '{}' and averaging method '{}' (mean value: {})".format(
                    observable_source_bands[observable_idx],
                    file_path,
                    observable_names[observable_idx],
                    observable_transforms[observable_idx],
                    observable_averaging_methods[observable_idx],
                    np.nanmean(source_data_interp)))
            # statistics over samples
            temp_transformed_sampled_data = transform_function(
                stats_on_all_samples(
                    source_data_interp,
                    pixel_mesh_east,
                    pixel_mesh_north,
                    sample_axis_east,
                    sample_size_east,
                    sample_axis_north,
                    sample_size_north,
                    additional_sampling_polygons,
                    method='nan_mean',
                ),
                observable_ranges[observable_idx],
                observable_transforms[observable_idx])
            # fill out the table
            observable_table[:,observable_idx] = np.kron(temp_transformed_sampled_data,np.ones(number_of_stacks))
            
    
    # mark rows in observable data table that have negative identifiers, nan-valued sar observables, infinite sar observables, or negative agb values
    invalid_rows = np.any(identifier_table<0,axis=1) | \
        np.any(np.isnan(observable_table[:,observable_is_required]),axis=1) | \
        np.any(~np.isfinite(observable_table[:,observable_is_required]),axis=1)
    # exclude invalid rows
    observable_table = observable_table[~invalid_rows,:]
    identifier_table = identifier_table[~invalid_rows,:]
    # number of rows in data table
    number_of_rows_in_observable_table = observable_table.shape[0]
            
            
    ### PREPARING PARAMETER TABLES
    parameter_property_names = ['lower_limit','upper_limit','initial_value']+['estimate_%d' % (ii) for ii in np.arange(number_of_subsets)+1]
    parameter_position_names = ['row_'+parameter_name for parameter_name in parameter_names]
    parameter_tables = []
    parameter_table_columns = []
    parameter_position_table = np.nan * np.zeros((
        number_of_rows_in_observable_table,number_of_parameters))
    # creating parameter matrices
    for parameter_idx,parameter_variability in enumerate(parameter_variabilities):
        # take out only the relevant identifiers (the ones that change as per parameter variability)
        # and create column names by adding four additional columns: min, max and initial value, and estimated value (later set to NaN)
        parameter_table_columns.append(np.concatenate((np.array(identifier_names)[parameter_variability],
                                                       np.array(parameter_property_names))))
        # create the minimal ID table (without unnecessary columns for those dimension across which the parameter doesn't change)
        temp_ids_table = identifier_table[:,np.where(parameter_variability)[0]]
        # create the last four columns
        temp_minmax_table = np.array([parameter_limits[parameter_idx]]) * np.ones((number_of_rows_in_observable_table,1))
        temp_initial_table = np.mean(temp_minmax_table,axis=1)
        temp_estimated_table = np.kron(np.ones((1,number_of_subsets)),np.array([np.mean(temp_minmax_table,axis=1)]).transpose())
        # create the full table
        # note: this table has initially the same shape as the observable table
        temp_full_table = np.column_stack((
                temp_ids_table,
                temp_minmax_table,
                temp_initial_table,
                temp_estimated_table))
        # take out unique rows and the inverse vector recreating the rows of the observable table
        # the inverse vector is critical as it relates the positions in the observable table to positions in each parameter table
        temp_full_table,temp_position_in_observable_table = np.unique(temp_full_table,axis=0,return_inverse=True)
        # set the last colum of the full table to nan (estimated value unavailable now)
        temp_full_table[:,-number_of_subsets:] = np.nan
        # append the table
        parameter_tables.append(temp_full_table)
        # include the inverse vector in the observable data table at the correct position
        parameter_position_table[:,parameter_idx] = temp_position_in_observable_table
                
    return (
        observable_table,
        observable_names,
        identifier_table,
        identifier_names,
        parameter_position_table,
        parameter_position_names,
        parameter_tables,
        parameter_table_columns,
        )