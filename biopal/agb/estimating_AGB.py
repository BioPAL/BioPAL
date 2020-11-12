# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 14:48:46 2020

@author: macie
"""

import numpy as np
import scipy as sp

from biopal.agb.processing_AGB import mean_on_rois


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
def stats_on_cal_polygons(data,pixel_axis_east,pixel_axis_north,cal_polygons,method):
    # here, create a function that calculates statistic in method on CAL data polygons
    #
    print('not implemented')
    return []
# %% function for calculating a given statistic on all samples on a grid and all polygons
def stats_on_all_samples(data,pixel_axis_east,pixel_axis_north,sampling_axis_east,sample_size_east,sampling_axis_north,sample_size_north,cal_polygons,method):
    stats = mean_on_rois(data,pixel_axis_east,pixel_axis_north,sampling_axis_east,sample_size_east,sampling_axis_north,sample_size_north,method)
    if cal_polygons:
        stats_cal_polygons = stats_on_cal_polygons(data,pixel_axis_east,pixel_axis_north,cal_polygons,method)
        stats = np.concatenate((stats,stats_cal_polygons))
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
         