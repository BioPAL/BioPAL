# -*- coding: utf-8 -*-
#%%
"""
Created on Tue Nov 10 14:48:46 2020

@author: macie
"""

import numpy as np
import scipy as sp
import logging
import os
import osr
import ogr
import gdal
from shapely.geometry import Polygon

from biopal.agb.processing_AGB import (
    check_intersection,
    interp2d_wrapper,
    merge_agb_intermediate,
    compute_processing_blocs_order,
)

# %%
def sample_and_tabulate_data(
    block_extents,  # extent of the current area for which the table is created
    pixel_axis_east,  # east north axes onto which data are interpolated
    pixel_axis_north,
    sampling_polygons,  # additional arbitrarily shaped polygons
    stack_in_block,  # flags whether each stack is in current block
    stack_info_table,  # info table with stack properties (stack id, headings, etc., defining the acquisition parameters)
    stack_info_table_column_names,  # column names for the abovementioned table
    formula_observables,
    observable_for_forest_class,
    formula_parameters,
    number_of_subsets,  # number of subsets to use (used to allocate columns in parameter tables)
):
    
    
    """
    (
        observable_table,
        observable_names
        identifier_table,
        identifier_names,
        parameter_position_table,
        parameter_position_names,
        parameter_tables,
        parameter_table_columns,
        sample_info_table,
        sample_info_table_columns,
    ) = sample_and_tabulate_data(
                                    block_extents,  # extent of the current area for which the table is created
                                    pixel_axis_east,  # east north axes onto which data are interpolated
                                    pixel_axis_north,
                                    sampling_polygons,  # additional arbitrarily shaped polygons
                                    stack_in_block,  # flags whether each stack is in current block
                                    stack_info_table,  # info table with stack properties (stack id, headings, etc., defining the acquisition parameters)
                                    stack_info_table_column_names,  # column names for the abovementioned table
                                    formula_observables,
                                    observable_for_forest_class,
                                    formula_parameters,
                                    number_of_subsets,  # number of subsets to use (used to allocate columns in parameter tables)
                                )
        
    Reads input data images, averages over polygons, and tabulates


    INPUT:
        block_extents contains the extents of the current block
        pixel_axis_east, pixel_axis_north are vectors defining the grid to which the read data has to be interpolated
        sampling_polygons is a list of Polygon objects over which the sampling is to be done
        stack_in_block is a vector of bools with the same length as the number of stacks, indicating whether each stack exists within the current block
        stack_info_table is a 2D array linking stack IDs to selected identifiers
        stack_info_table_column_names is the corresponding list of column names
        formula_observables is a named tuple with information about formula observables
        observale_for_forest_class is the name of the observable that acts as a forest_class observable
        formula_parameters is a named tuple with information about formula parameters
        number_of_subsets is the number subsets (independent tests) to be run
        
    OUTPUT:
        observable_table is a 2D table with observable data averaged over the polygons
        observable_names is a corresponding list of observables (column names for the table above)
        identifier_table is a 2D table with identifiers matching the observable_table
        identifier_names is a corresponding list of identifier names (column names for the table above)
        parameter_position_table is a 2D table with columnwise unique indices for each parameter
        parameter_position_names is the corresponding list of column names for the table above
        parameter_tables is a list of minimal-sized 2D arrays within which the estimated parameter values for each parameter will be saved (one list element for each unknown parameter)
        parameter_table_columns is a corresponding list of lists with table column names
        sample_info_table is a 2D table with information about each sample
        sample_info_table_columns is the corresponding list with names for the columns of the table above
        
        
        
        
        
    """

    # derived parameters
    number_of_stacks = stack_info_table.shape[0]
    stack_id_vector = stack_info_table[:, 0]
    number_of_samples = len(sampling_polygons)
    number_of_parameters = len(formula_parameters.name)
    number_of_observables = len(formula_observables.name)
    pixel_mesh_east, pixel_mesh_north = np.meshgrid(pixel_axis_east, pixel_axis_north)


    # defining names for identifiers (sampleID & forest class ID and then all columns from stack info table)
    identifier_names = ["sample_id", "forest_class_id"] + stack_info_table_column_names #+ ["resolution_m"]
    number_of_identifiers = len(identifier_names)
    identifier_table = np.nan * np.zeros((number_of_samples * number_of_stacks, number_of_identifiers))

    # filling out the first column with sample IDs
    identifier_table[:, 0] = np.kron(np.arange(number_of_samples), np.ones(number_of_stacks))

    # filling out columns 3-8 with stack IDs and corresponding other identifications from stack_info_table
    identifier_table[:, 2] = np.kron(np.ones(number_of_samples), stack_id_vector)
    for id_idx in range(5):
        identifier_table[:, 3 + id_idx] = sp.interpolate.interp1d(
            stack_id_vector, stack_info_table[:, 1 + id_idx], kind="nearest"
        )(identifier_table[:, 2])
        
    
    # allocate observable table
    observable_table = np.nan * np.zeros((number_of_samples * number_of_stacks , number_of_observables))
    sample_info_table_columns = ['easting','northing','area','east_range','north_range']
    sample_info_table = np.nan*np.zeros((number_of_samples*number_of_stacks,len(sample_info_table_columns)))
    
    
    ## get info about samples (coordinates, area, shape factor)
    # calculating sample area and coordinates
    sample_info_table = []
    for current_map,current_stat,current_name in zip(
            [pixel_mesh_east,
             pixel_mesh_north,
             pixel_mesh_east*0 + np.abs(np.diff(pixel_axis_east)[0]*np.diff(pixel_axis_north)[0]),
             pixel_mesh_east,
             pixel_mesh_north],
            ["mean","mean","sum","range","range"],
            sample_info_table_columns):
        
        logging.info("AGB: calculating sample statistic '{}'".format(current_name))
        
        sample_info_table.append(
            np.round(
                stats_on_all_samples(
                    current_map,
                    0,
                    pixel_mesh_east,
                    pixel_mesh_north,
                    sampling_polygons,
                    current_stat,
                )
            )
        )
    
    
    sample_info_table = np.column_stack(sample_info_table)
    sample_info_table = np.kron(sample_info_table, np.ones((number_of_stacks,1)))
    
    
    ### READING OBSERVABLE DATA
    # cycle through observable sources in the stack
    for observable_idx in range(number_of_observables):
        # number of stacks in current observable (must be either 1 or number_of_stacks)
        current_number_of_stacks = len(formula_observables.source_paths[observable_idx])

        is_stacked = (current_number_of_stacks==number_of_stacks)

        # cycle over stacks
        for stack_idx in range(current_number_of_stacks):

            # go ahead only if current stack is (at least partially) contained in the current parameter block:
            if (is_stacked and stack_in_block[stack_idx]) or not is_stacked:

                # cycle over all the equi7 tiles, interpolate over pixel grid and average
                source_data_interp = np.NaN * np.zeros((len(pixel_axis_north), len(pixel_axis_east)), dtype="float")

                for file_idx, file_source in enumerate(formula_observables.source_paths[observable_idx][stack_idx]):

                    if file_source[0]!='none':#np.round(file_source[2])<=current_resolution:
                        source_data_interp_curr = interp2d_wrapper(
                            file_source[0], file_source[1] + 1, pixel_axis_east, pixel_axis_north, fill_value=np.NaN,
                        )

                        source_data_interp = merge_agb_intermediate(
                            source_data_interp, source_data_interp_curr, method="nan_mean"
                        )

                        # # masking the stack:
                        # source_data_interp[forest_class_map_interp == 0] = np.NaN
                        if is_stacked:
                            logging.info(
                                "AGB: sampling and transforming data for stacked observable '{}' (stack {}/{}, file {}/{})".format(
                                    formula_observables.name[observable_idx],
                                    stack_idx + 1,
                                    number_of_stacks,
                                    file_idx + 1,
                                    len(formula_observables.source_paths[observable_idx][stack_idx]),
                                )
                            )
                        else:
                            logging.info(
                                "AGB: sampling and transforming data for unstacked observable '{}' (file {}/{})".format(
                                    formula_observables.name[observable_idx],
                                    file_idx + 1,
                                    len(formula_observables.source_paths[observable_idx][stack_idx]),
                                )
                            )

                    # calculate sample statistics
                    temp_transformed_sampled_data = transform_function(
                        stats_on_all_samples(
                            source_data_interp,
                            formula_observables.source_resolution[observable_idx],
                            pixel_mesh_east,
                            pixel_mesh_north,
                            sampling_polygons,
                            formula_observables.averaging_method[observable_idx],
                        ),
                        formula_observables.limits[observable_idx],
                        formula_observables.transform[observable_idx],
                    )
                    
                    # check if observable is stacked (i.e., the list length equals number of stacks)
                    if is_stacked:
                        # find rows where the stack id matches the current stack id
                        current_rows = (identifier_table[:, 2] == stack_id_vector[stack_idx])# & (identifier_table[:,-1]==current_resolution)
                        # fill out the table
                        observable_table[current_rows, observable_idx] = temp_transformed_sampled_data
                    else:
                        
                        # fill out the table by replicating the averaged data
                        observable_table[:, observable_idx] = np.kron(temp_transformed_sampled_data, np.ones(number_of_stacks))
                        
                        
        # break if this is a required observable and all sampled data are nan
        # (speeds up the reading)
        if np.all(np.isnan(observable_table[:, observable_idx])) & formula_observables.is_required[observable_idx]:
            logging.info("AGB: no data for the first required observable, skipping the rest")
            break

    
    # repeating it across stacks and inserting into observable data table (note: forest class assumed constant across stacks)
    # identifier_table[:, 1] = np.kron(temp_forest_class_vector, np.ones(number_of_stacks))
    forest_class_position = np.where(match_string_lists(formula_observables.name,[observable_for_forest_class]).flatten()>=0)[0]
    if len(forest_class_position)==0:
        logging.error('AGB: cannot find forest class observable {}.' .format(observable_for_forest_class))
    else:
        forest_class_position = forest_class_position[0]
    identifier_table[:,1] = observable_table[:,forest_class_position]
    

    # observable_table = np.nan*observable_table # this is just a dummy thing to check the behaviour of this function in case of lack of data

    # mark rows in observable data table that have negative identifiers, nan-valued sar observables, infinite sar observables, or negative agb values
    invalid_rows = (
        np.any(identifier_table < 0, axis=1)
        | np.any(np.isnan(observable_table[:, formula_observables.is_required]), axis=1)
        | np.any(~np.isfinite(observable_table[:, formula_observables.is_required]), axis=1)
    )
    # exclude invalid rows
    observable_table = observable_table[~invalid_rows, :]
    identifier_table = identifier_table[~invalid_rows, :]
    sample_info_table = sample_info_table[~invalid_rows, :]
    # number of rows in data table
    number_of_rows_in_observable_table = observable_table.shape[0]

    ### PREPARING PARAMETER TABLES
    parameter_property_names = ["lower_limit", "upper_limit", "initial_value"] + [
        "estimate_%d" % (ii) for ii in np.arange(number_of_subsets) + 1
    ]
    parameter_position_names = ["row_" + parameter_name for parameter_name in formula_parameters.name]
    parameter_tables = []
    parameter_table_columns = []
    parameter_position_table = np.nan * np.zeros((number_of_rows_in_observable_table, number_of_parameters))
    # creating parameter matrices
    for parameter_idx, parameter_variability in enumerate(formula_parameters.parameter_variabilities):
        # take out only the relevant identifiers (the ones that change as per parameter variability)
        # and create column names by adding four additional columns: min, max and initial value, and estimated value (later set to NaN)
        parameter_table_columns.append(
            np.concatenate((np.array(identifier_names)[parameter_variability], np.array(parameter_property_names)))
        )
        # create the minimal ID table (without unnecessary columns for those dimension across which the parameter doesn't change)
        temp_ids_table = identifier_table[:, np.where(parameter_variability)[0]]
        # create the last four columns
        temp_minmax_table = np.array([formula_parameters.limits[parameter_idx]]) * np.ones(
            (number_of_rows_in_observable_table, 1)
        )
        temp_initial_table = np.mean(temp_minmax_table, axis=1)
        temp_estimated_table = np.kron(
            np.ones((1, number_of_subsets)), np.array([np.mean(temp_minmax_table, axis=1)]).transpose()
        )
        # create the full table
        # note: this table has initially the same shape as the observable table
        temp_full_table = np.column_stack((temp_ids_table, temp_minmax_table, temp_initial_table, temp_estimated_table))
        # take out unique rows and the inverse vector recreating the rows of the observable table
        # the inverse vector is critical as it relates the positions in the observable table to positions in each parameter table
        temp_full_table, temp_position_in_observable_table = np.unique(temp_full_table, axis=0, return_inverse=True)
        # set the last colum of the full table to nan (estimated value unavailable now)
        temp_full_table[:, -number_of_subsets:] = np.nan
        # append the table
        parameter_tables.append(temp_full_table)
        # include the inverse vector in the observable data table at the correct position
        parameter_position_table[:, parameter_idx] = temp_position_in_observable_table

    return (
        observable_table,
        formula_observables.name,
        identifier_table,
        identifier_names,
        parameter_position_table,
        parameter_position_names,
        parameter_tables,
        parameter_table_columns,
        sample_info_table,
        sample_info_table_columns,
    )


# %%


def fit_formula_to_random_subsets(
    formula_terms,
    number_of_subsets,
    observable_table,
    observable_names,
    identifier_table,
    identifier_names,
    parameter_position_table,
    parameter_names,
    parameter_tables,
    parameter_table_columns,
    parameter_variabilities,
    calibration_fraction,
    estimation_fraction,
    calibration_areas_per_test,  # proc_conf.AGB.min_number_of_cals_per_test
    estimation_areas_per_test,  # proc_conf.AGB.min_number_of_rois_per_test
    transfer_function_name,
):

    
    """
    (
        parameter_tables,
        space_invariant_parameter_table,
        space_invariant_parameter_names,
        space_variant_parameter_table,
        space_variant_parameter_names,
    ) = fit_formula_to_random_subsets(
            formula_terms,
            number_of_subsets,
            observable_table,
            observable_names,
            identifier_table,
            identifier_names,
            parameter_position_table,
            parameter_names,
            parameter_tables,
            parameter_table_columns,
            parameter_variabilities,
            calibration_fraction,
            estimation_fraction,
            calibration_areas_per_test, 
            estimation_areas_per_test, 
            transfer_function_name,
        )
    
    
    Fits space invariant and space variant parameters to the data by minimising the residuals defined in the formulas
    
    
    INPUT: 
        formula_terms is a named tuple with formula strings and associated formula weights
        number_of_subsets is the number of random subsets to be created
        observable_table is a 2D numpy array with observable data
        observable_names is a corresponding list of observable names
        identifier_table is a 2D numpy array with identifiers
        identifier_names is a corresponding list of identifier names
        parameter_position_table is a 2D numpy array with indices indicating how the parameters to be estimated vary across samples
        parameter_names is a corresonding list of parameter names
        parameter_tables is a list of parameter tables, each containing the minimal number of rows and columns for each parameter, where the rows contain individual parameter values and columns contain estimates from each subset and required identifiers
        parameter_columns is a list of lists with parameter column names
        parameter_variabilities contains the variabilities of each parameter across the predefined dimensions
        calibration_fraction is the fraction of samples with all required and non-requried data to be used in each subset
        estimation_fraction is the fraction of all samples with all required, but with some missing non-required data to be used in each subset
        calibration_areas_per_test is the minimal number of calibration areas required for this function to be run
        estimation_areas_per_test is the minimal number of esti9mation areas required for this function to be run
        transfer_function_name is the name of the transfer function used in this function
        
    OUTPUT:
        parameter_tables is the updated list of parameter tables
        space_invariant_parameter_table is a 2D numpy array on the same format as observable_table (same number of rows), with the estimated space invariant parameter values replicated as parameter_position_table
        space_invariant_parameter_names is a corresponding list of space invariant parameter names
        space_variant_parameter_table is a 2D numpy array on the same format as observable_table (same number of rows), with the estimated space variant parameter values replicated as parameter_position_table
        space_variant_parameter_names is a corresponding list of space variant parameter names
 
    
    
    """
    
    
    ### CREATE CALIBRATION AND ESTIMATION SUBSETS
    logging.info("AGB: creating {} subsets".format(number_of_subsets))

    # select rows with available agb information as calibration data and those without as estimation data
    columns_with_nans = np.all(np.isnan(observable_table),axis=0)
    calibration_rows = np.where(np.all(~np.isnan(observable_table[:,~columns_with_nans]), axis=1))[0]
    estimation_rows = np.where(np.any(np.isnan(observable_table[:,~columns_with_nans]), axis=1))[0]
    calibration_sample_ids = np.unique(identifier_table[calibration_rows, 0])
    estimation_sample_ids = np.unique(identifier_table[estimation_rows, 0])

    # calculate subset sizes
    estimation_subset_size = np.int32(np.ceil(len(estimation_sample_ids) * estimation_fraction))
    calibration_subset_size = np.int32(np.ceil(len(calibration_sample_ids) * calibration_fraction))

    # find random data subsetting vectors making sure that the number of calibration and estimation areas
    # is the same in all
    subset_indexing_vectors = []
    number_of_accepted_subsets = 0
    max_number_of_subsets_to_test = number_of_subsets * 10
    tested_subset_counter = 0
    while (number_of_accepted_subsets < number_of_subsets) & (tested_subset_counter < max_number_of_subsets_to_test):

        # create a random subset of calibration and estimation samples
        current_random_estimation_subset = np.sort(
            np.random.permutation(estimation_sample_ids)[:estimation_subset_size]
        )
        current_random_calibration_subset = np.sort(
            np.random.permutation(calibration_sample_ids)[:calibration_subset_size]
        )

        # calculate the minimal number of calibration and estimation samples for the space-invariant parameters
        # (for the latter, we use the column with parameter positions in parameter tables - the same value indicates the same parameter)
        current_parameter_position_columns = np.where(~np.row_stack(parameter_variabilities)[:, 0])[0]

        current_calibration_rows = np.isin(identifier_table[:, 0], current_random_calibration_subset)
        min_number_of_calibration_measurements_per_space_invariant_parameter = np.inf
        for column_idx in current_parameter_position_columns:
            # calculate the minimal number of samples for all parameter values within this column and all different parameters until the current one
            if np.any(current_calibration_rows):
                min_number_of_calibration_measurements_per_space_invariant_parameter = np.minimum(
                    min_number_of_calibration_measurements_per_space_invariant_parameter,
                    np.min(
                        np.unique(parameter_position_table[current_calibration_rows, column_idx], return_counts=True)[1]
                    ),
                )
            else:
                min_number_of_calibration_measurements_per_space_invariant_parameter = 0

        current_estimation_rows = np.isin(identifier_table[:, 0], current_random_estimation_subset)
        min_number_of_estimation_measurements_per_space_invariant_parameter = np.inf
        for column_idx in current_parameter_position_columns:
            # calculate the minimal number of samples for all parameter values within this column and all different parameters until the current one
            if np.any(current_estimation_rows):
                min_number_of_estimation_measurements_per_space_invariant_parameter = np.minimum(
                    min_number_of_estimation_measurements_per_space_invariant_parameter,
                    np.min(
                        np.unique(parameter_position_table[current_estimation_rows, column_idx], return_counts=True)[1]
                    ),
                )
            else:
                min_number_of_estimation_measurements_per_space_invariant_parameter = 0

        # if the minimal number of samples is larger than the one specified in the xml configuration file, accept this subset
        # (at the moment, we don't perform other tests, which means that subsets may be repeated)
        if (min_number_of_calibration_measurements_per_space_invariant_parameter >= calibration_areas_per_test) & (
            min_number_of_estimation_measurements_per_space_invariant_parameter >= estimation_areas_per_test
        ):
            subset_indexing_vectors.append(
                np.isin(
                    identifier_table[:, 0],
                    np.sort(np.concatenate((current_random_calibration_subset, current_random_estimation_subset))),
                )
            )
            number_of_accepted_subsets += 1

        tested_subset_counter += 1
    if number_of_accepted_subsets == 0:
        logging.info("AGB: no accepted subsets found in current block.")
        space_invariant_parameter_table = []
        space_invariant_parameter_names = []
        space_variant_parameter_table = []
        space_variant_parameter_names = []

    else:
        if number_of_accepted_subsets < number_of_subsets:
            logging.warning(
                "AGB: number of accepted subsets ({}) is less than the required number of subsets ({}).".format(
                    number_of_accepted_subsets, number_of_subsets
                )
            )
            
            
        # detect observables without data (used to remove formula terms that would generate nans)
        observables_without_data = [observable_name for is_nan,observable_name in zip(np.all(np.isnan(observable_table),axis=0),observable_names) if is_nan]
        terms_with_nan_observables = np.any(match_string_lists(formula_terms.string,observables_without_data)>=0,axis=1)
        
        if np.any(terms_with_nan_observables):
            logging.warning(
                "AGB: skipping formula terms: {} in steps 1-3 due to lack of useful data for observables: {}.".format(
                    ', '.join(['%d' % (ii+1) for ii in np.where(terms_with_nan_observables)[0]]),
                    ', '.join(observables_without_data)))
        
        
        ### FIRST, WE PERFORM THE ESTIMATION OF PARAMETERS FOR SUBSETS
        
        
        logging.info("AGB: parameter estimation step 1: estimation for subsets")
        
        
        #### select relevant formula and weights for this step
        terms_with_zero_weight_fitting = np.array(formula_terms.formula_weights.fitting) == 0
        if np.any(terms_with_zero_weight_fitting):
            logging.warning(
                "AGB: skipping formula terms: {} in step 1 due to zero weights.".format(
                    ', '.join(['%d' % (ii+1) for ii in np.where(terms_with_zero_weight_fitting)[0]])))
        
        terms_to_take_fitting = ~(terms_with_nan_observables | terms_with_zero_weight_fitting)
        formula_fitting = subset_iterable(formula_terms.string,terms_to_take_fitting)
        formula_weights_fitting = subset_iterable(formula_terms.formula_weights.fitting,terms_to_take_fitting)
        
        # find observables and parameters in formula and create
        # vectors for selecting parameters and observables that exist in formula
        observables_in_formula_fitting = np.any(match_string_lists(formula_fitting, observable_names) >= 0, axis=0)
        observables_in_parameters = np.any(match_string_lists(parameter_names, observable_names) >= 0, axis=0)
        parameters_in_formula_fitting = np.any(match_string_lists(formula_fitting, parameter_names) >= 0, axis=0)
        # find parameters that do not change between samples
        space_invariant_parameters = False == np.column_stack(parameter_variabilities)[0, :]
        # some more flag vectors
        space_invariant_parameters_in_formula_fitting = space_invariant_parameters & parameters_in_formula_fitting
        
        # loop through calibration subsets
        for subset_idx, current_subset in enumerate(subset_indexing_vectors):

            logging.info(
                "AGB: minimising cost function for subset {} out of {}...".format(subset_idx + 1, number_of_subsets)
            )

            # subset parameter and observable tables
            current_parameter_position_table = parameter_position_table[current_subset, :][:, parameters_in_formula_fitting]
            current_parameter_names = subset_iterable(parameter_names, parameters_in_formula_fitting, False)
            current_observable_table = observable_table[current_subset, :][:, observables_in_formula_fitting]
            current_observable_names = subset_iterable(observable_names, observables_in_formula_fitting, False)
            # create min-max tables
            individual_parameter_min_max_tables = subset_iterable(parameter_tables, parameters_in_formula_fitting, False)
            for parameter_idx in range(len(individual_parameter_min_max_tables)):
                individual_parameter_min_max_tables[parameter_idx] = individual_parameter_min_max_tables[parameter_idx][
                    :, -number_of_subsets - 3 : -number_of_subsets - 1
                ]

            # estimate both parameters and AGB for the subset
            (current_lut_all_parameters, _, _,) = fit_formula_to_table_data(
                formula_fitting,
                formula_weights_fitting,
                current_observable_table,#subset_iterable(current_observable_table.transpose(),~np.all(np.isnan(current_observable_table),axis=0),return_array=True).transpose(),
                current_observable_names,#subset_iterable(current_observable_names,~np.all(np.isnan(current_observable_table),axis=0)),
                current_parameter_position_table,
                current_parameter_names,
                individual_parameter_min_max_tables,
                transfer_function_name,
            )

            # fill out parameter tables with estimates of space invariant parameters
            for current_parameter_idx in np.where(space_invariant_parameters_in_formula_fitting)[0]:
                # identify column in the current output table
                current_column_idx = np.where(
                    np.array(current_parameter_names) == np.array(parameter_names[current_parameter_idx])
                )[0][0]
                # identify valid rows in the current output lut
                current_rows = (current_lut_all_parameters[:, 1] == current_column_idx) & (
                    np.abs(current_lut_all_parameters[:, -2] - current_lut_all_parameters[:, -1])
                    > 1e-10  # checking if the parameter has changed from the initial value
                )
                # write the relevant rows of the parameter table
                parameter_tables[current_parameter_idx][
                    np.int32(current_lut_all_parameters[current_rows, 2]), -number_of_subsets + subset_idx
                ] = current_lut_all_parameters[current_rows, -1]






        ### THEN, WE ESTIMATE SPACE VARIANT PARAMETERS FOR ALL SAMPLES USING SPACE INVARIANT PARAMETERS FROM SUBSETS
        logging.info("AGB: parameter estimation step 2: estimation of AGB from subset estimated parameters")
        
        
        #### select relevant formula and weights for this step
        terms_with_zero_weight_estimation1 = np.array(formula_terms.formula_weights.estimation1) == 0
        if np.any(terms_with_zero_weight_estimation1):
            logging.warning(
                "AGB: skipping formula terms: {} in step 2 due to zero weights.".format(
                    ', '.join(['%d' % (ii+1) for ii in np.where(terms_with_zero_weight_estimation1)[0]])))
        terms_to_take_estimation1 = ~(terms_with_nan_observables | terms_with_zero_weight_estimation1)
        formula_estimation1 = subset_iterable(formula_terms.string,terms_to_take_estimation1)
        formula_weights_estimation1 = subset_iterable(formula_terms.formula_weights.estimation1,terms_to_take_estimation1)
        
        observables_in_formula_estimation1 = np.any(match_string_lists(formula_estimation1, observable_names) >= 0, axis=0)
        parameters_in_formula_estimation1 = np.any(match_string_lists(formula_estimation1, parameter_names) >= 0, axis=0)
        space_invariant_parameters_in_formula_estimation1 = space_invariant_parameters & parameters_in_formula_estimation1
        space_variant_parameters_in_formula_estimation1 = ~space_invariant_parameters & parameters_in_formula_estimation1
        observables_in_formula_estimation1_not_in_parameters = observables_in_formula_estimation1 & ~observables_in_parameters

        
        # loop through calibration subsets
        for subset_idx in range(number_of_subsets):

            logging.info(
                "AGB: estimating space-variant parameters for all samples using space-invariant parameter set {} out of {}...".format(
                    subset_idx + 1, number_of_subsets
                )
            )

            # create a table with all space invariant parameters
            space_invariant_parameter_names = subset_iterable(
                parameter_names, space_invariant_parameters_in_formula_estimation1, False
            )
            current_space_invariant_parameter_table = []
            for current_parameter_idx in np.where(space_invariant_parameters_in_formula_estimation1)[0]:
                # extract current parameter table
                current_parameter_table = parameter_tables[current_parameter_idx]
                # extract the current position vector
                current_parameter_position_vector = np.int32(parameter_position_table[:, current_parameter_idx])
                # extract the column with the current subset estimates
                current_column_in_parameter_table = -number_of_subsets + subset_idx
                # add the current data
                current_space_invariant_parameter_table.append(
                    current_parameter_table[current_parameter_position_vector, current_column_in_parameter_table]
                )
            current_space_invariant_parameter_table = np.column_stack(current_space_invariant_parameter_table)

            # these parameters are now treated as observables, so
            # space invariant parameter table is merged with observable table
            new_observable_table = np.column_stack(
                (observable_table[:, observables_in_formula_estimation1_not_in_parameters], current_space_invariant_parameter_table)
            )
            new_observable_names = (
                subset_iterable(observable_names, observables_in_formula_estimation1_not_in_parameters, False)
                + space_invariant_parameter_names
            )

            # now, only the space-variant parameters are treated us unknown parameters
            # the corresponding tables and lists are now created
            new_parameter_position_table = parameter_position_table[:, space_variant_parameters_in_formula_estimation1]
            new_parameter_names = subset_iterable(parameter_names, space_variant_parameters_in_formula_estimation1, False)
            new_individual_parameter_min_max_tables = subset_iterable(
                parameter_tables, space_variant_parameters_in_formula_estimation1, False
            )
            for parameter_idx in range(len(new_individual_parameter_min_max_tables)):
                new_individual_parameter_min_max_tables[parameter_idx] = new_individual_parameter_min_max_tables[
                    parameter_idx
                ][:, -number_of_subsets - 3 : -number_of_subsets - 1]

            # estimate space variant parameters for all samples
            (current_lut_space_variant_parameters, space_variant_parameter_table, _,) = fit_formula_to_table_data(
                formula_estimation1,
                formula_weights_estimation1,
                new_observable_table,
                new_observable_names,
                new_parameter_position_table,
                new_parameter_names,
                new_individual_parameter_min_max_tables,
                transfer_function_name,
            )

            # fill out parameter tables with estimates of space invariant parameters
            for current_parameter_idx in np.where(space_variant_parameters_in_formula_estimation1)[0]:
                # identify column
                current_column_idx = np.where(
                    np.array(new_parameter_names) == np.array(parameter_names[current_parameter_idx])
                )[0][0]
                # identify valid rows
                current_rows = (current_lut_space_variant_parameters[:, 1] == current_column_idx) & (
                    np.abs(current_lut_space_variant_parameters[:, -2] - current_lut_space_variant_parameters[:, -1])
                    > 1e-10
                )
                # save current estimates to parameter tables
                parameter_tables[np.int32(current_parameter_idx)][
                    np.int32(current_lut_space_variant_parameters[current_rows, 2]), -number_of_subsets + subset_idx
                ] = current_lut_space_variant_parameters[current_rows, -1]

        ### FINALLY, WE ESTIMATE SPACE INVARIANT PARAMETERS USING THE AVERAGE SPACE VARIANT PARAMETER VALUES FROM ALL SUBSETS FOR ALL SAMPLES

        logging.info("AGB: parameter estimation step 3: fitting space-invariant parameters using space-variant parameter estimate")

        terms_with_zero_weight_estimation2 = np.array(formula_terms.formula_weights.estimation2) == 0
        if np.any(terms_with_zero_weight_estimation2):
            logging.warning(
                "AGB: skipping formula terms: {} in step 3 due to zero weights.".format(
                    ', '.join(['%d' % (ii+1) for ii in np.where(terms_with_zero_weight_estimation2)[0]])))
        terms_to_take_estimation2 = ~(terms_with_nan_observables | terms_with_zero_weight_estimation2)
        formula_estimation2 = subset_iterable(formula_terms.string,terms_to_take_estimation2)
        formula_weights_estimation2 = subset_iterable(formula_terms.formula_weights.estimation2,terms_to_take_estimation2)
        
        observables_in_formula_estimation2 = np.any(match_string_lists(formula_estimation2, observable_names) >= 0, axis=0)
        parameters_in_formula_estimation2 = np.any(match_string_lists(formula_estimation2, parameter_names) >= 0, axis=0)
        space_invariant_parameters_in_formula_estimation2 = space_invariant_parameters & parameters_in_formula_estimation2
        space_variant_parameters_in_formula_estimation2 = ~space_invariant_parameters & parameters_in_formula_estimation2
        observables_in_formula_estimation2_not_in_parameters = observables_in_formula_estimation2 & ~observables_in_parameters


        # names for table columns
        space_variant_parameter_names = subset_iterable(parameter_names, space_variant_parameters_in_formula_estimation2, False)
        # create table with all space invariant parameters
        current_space_variant_parameter_table = []
        for position_in_parameter_table_list in np.where(space_variant_parameters_in_formula_estimation2)[0]:
            current_parameter_table = parameter_tables[position_in_parameter_table_list]
            current_parameter_position_vector = np.int32(parameter_position_table[:, position_in_parameter_table_list])
            current_columns_in_parameter_table = np.arange(-number_of_subsets, 0)
            current_space_variant_parameter_table.append(
                np.mean(
                    current_parameter_table[current_parameter_position_vector, :][
                        :, current_columns_in_parameter_table
                    ],
                    axis=1,
                )
            )
        current_space_variant_parameter_table = np.column_stack(current_space_variant_parameter_table)

        # new observable table is the combination of observable table and space variant parameter table
        new_observable_table = np.column_stack(
            (observable_table[:, observables_in_formula_estimation2_not_in_parameters], current_space_variant_parameter_table)
        )
        new_observable_names = (
            subset_iterable(observable_names, observables_in_formula_estimation2_not_in_parameters, False)
            + space_variant_parameter_names
        )

        new_parameter_position_table = parameter_position_table[:, space_invariant_parameters_in_formula_estimation2]
        new_parameter_names = subset_iterable(parameter_names, space_invariant_parameters_in_formula_estimation2, False)
        new_individual_parameter_min_max_tables = subset_iterable(
            parameter_tables, space_invariant_parameters_in_formula_estimation2, False
        )
        for parameter_idx in range(len(new_individual_parameter_min_max_tables)):
            new_individual_parameter_min_max_tables[parameter_idx] = new_individual_parameter_min_max_tables[
                parameter_idx
            ][:, -number_of_subsets - 3 : -number_of_subsets - 1]

        # estimate space variant parameters for all samples
        (_, space_invariant_parameter_table, _,) = fit_formula_to_table_data(
            formula_estimation2,
            formula_weights_estimation2,
            new_observable_table,
            new_observable_names,
            new_parameter_position_table,
            new_parameter_names,
            new_individual_parameter_min_max_tables,
            transfer_function_name,
        )
    return (
        parameter_tables,
        space_invariant_parameter_table,
        space_invariant_parameter_names,
        space_variant_parameter_table,
        space_variant_parameter_names,
    )


# %% functions needed for function above
# swap variable names in formulas to slices of an array
def swap_names_and_merge_formula(
    original_formulas, observable_names, parameter_names, new_table_name, use_observable_if_repeated_and_available=True
):
    
    """
    
    new_formulas = swap_names_and_merge_formula(
                                    original_formulas, observable_names, parameter_names, new_table_name, use_observable_if_repeated_and_available=True
                                )
    
    Swaps observable names and parameter names in the original formulas for slices of a table. Column names are determined by position of parameter/observable names in the combined list
    of observable and parameter names
    
    
    
    INPUT:
        original_formulas is a list with formula strings
        observable_names is a list of observable names
        parameter_names is a list of parameter names
        new_table_name is the name for the new table
        use_observable_if_repeated_and_available indicates whether observables should be prioritised if the same name occurs both in the parameter_names and observable_names lists
    
    OUTPUT:
        new_formulas is a list of formula strings where the observable/parameter names have been replaced with slices of the new table
    
    """
    original_variable_names = observable_names + parameter_names
    unique_variable_names, name_counts = np.unique(np.array(original_variable_names), return_counts=True)
    new_formula = []
    for current_formula in original_formulas:
        for unique_variable_name, name_count in zip(unique_variable_names, name_counts):
            if name_count == 1:
                position_in_variable_names_vector = np.where(np.array(original_variable_names) == unique_variable_name)[
                    0
                ][0]
            elif name_count == 2:
                position_in_variable_names_vector = np.where(np.array(original_variable_names) == unique_variable_name)[
                    0
                ][np.int32(~use_observable_if_repeated_and_available)]
            current_formula = current_formula.replace(
                unique_variable_name, new_table_name + ("[:,%d]" % (position_in_variable_names_vector))
            )
        new_formula.append(current_formula)
    return new_formula

# %%
# function for converting columnwise indices (which are repeated within the same column if they represent identical values,
# but which may be repeated across different columns without meaning that they represent identical values)
# to unique indices (which are only repeated within the same table if they are to have identical values)
def regularise_indices(columnwise_index_table):
    
    """
    
    regularised_index_table = regularise_indices(columnwise_index_table)
    
    Converts columnwise indices to global indices. 
    
    INPUT:
        columnwise_index_table is a 2D numpy array with different numbers indicating different parameter values within each column
        
    OUTPUT:
        regularised_index_table is a 2D numpy array with different numbers indicating different parameter values within the entire table
        
    
    
    """

    offset = 0
    unique_index_table = []  # position in a single x-vector
    columnwise_to_unique_index_lut = (
        []
    )  # lut for converting between parameter id and column and parameter position in x
    for column_idx, parameter_column in enumerate(columnwise_index_table.transpose()):
        # add -1 at the beginning to avoid vectors with one element (which will not work with interp1d)
        old_indices = np.concatenate((-1 * np.ones(1), np.unique(parameter_column)))
        # new indices is a simple sequence from 0 to number of parameters-1 + offset due to previous parameters
        new_indices = np.arange(len(old_indices)) - 1 + offset
        # convert parameter indices and add to the list
        unique_index_table.append(sp.interpolate.interp1d(old_indices, new_indices, kind="nearest")(parameter_column))
        # save the lut, removing the first, unnecessary element
        columnwise_to_unique_index_lut.append(
            np.column_stack(
                (
                    new_indices[old_indices > -1],
                    column_idx + np.zeros(len(old_indices[old_indices > -1])),
                    old_indices[old_indices > -1],
                )
            )
        )
        # update offset based on current parameter column
        offset = np.max(new_indices) + 1
    # convert the list of vectors to an array
    unique_index_table = np.int32(np.column_stack(unique_index_table))
    # stack all luts to one lut
    columnwise_to_unique_index_lut = np.row_stack(columnwise_to_unique_index_lut)
    # # length of beta vector
    # length_of_p_vector = columnwise_to_unique_index_lut.shape[0]
    return unique_index_table, columnwise_to_unique_index_lut

# %%
def cost_function(
    x_vector, converted_formulas, formula_weights, observable_table, index_table, name_of_table_in_converted_formula, transfer_function, return_one_value=True
):
    """
    cost_function_value = cost_function(
                                x_vector, converted_formulas, formula_weights, observable_table, index_table, name_of_table_in_converted_formula, transfer_function, return_one_value=True
                            )
    
    
    
    Calculates the cost function value from a one-dimensional vector
    
    
    INPUT:
        
        x_vector is a one-dimensional vector as required by the sp.optimize.minimize function
        converted_formulas is a list of formula strings converted to use columns of one single table
        formula_weights is a list of formula weights
        observable_table is a 2D numpy array with observables
        index_table is a 2D numpy array with parameter indices in x_vector
        name_of_table_in_converted_formula is a string with the name that is used in the converted formulas for the combined table
        transfer_function is a string with the transfer function type
        return_one_value is a bool indicating if one value should be returned (default: True), or if the function should return one value for each formula term (False)
        
    OUTPUT:
        
        cost_function_value  is either a scalar or an array of the same length as the list of formula strings
    
    
    
    
    """
    
    
    
    
    p_vector = transfer_function(x_vector)
    table_in_converted_formula = np.column_stack((observable_table, p_vector[index_table]))
    final_expression = '0'
    for converted_formula,formula_weight in zip(converted_formulas,formula_weights):
        final_expression += ',%.18f*np.nanmean((%s)**2)' % (formula_weight/np.sum(formula_weights)*len(formula_weights),converted_formula)
    final_expression = 'np.array([%s])' % (final_expression)
    total_cost = eval(final_expression, {name_of_table_in_converted_formula: table_in_converted_formula, "np": np})[1:]
    if return_one_value:
        return np.sqrt(np.nansum(total_cost))
    else:
        return total_cost


# %%

def fit_formula_to_table_data(
    original_formula,
    formula_weights,
    observable_table,
    observable_table_column_names,
    parameter_position_table,
    parameter_position_table_column_names,
    individual_parameter_min_max_tables,
    transfer_function_name,
):
   

    """
    
    
    estimated_parameter_lut, estimated_parameter_table, cost_function_values = fit_formula_to_table_data(
                                                                                            original_formula,
                                                                                            formula_weights,
                                                                                            observable_table,
                                                                                            observable_table_column_names,
                                                                                            parameter_position_table,
                                                                                            parameter_position_table_column_names,
                                                                                            individual_parameter_min_max_tables,
                                                                                            transfer_function_name,
                                                                                        )
    
    Performs minimisation of the cost function defined by a formula and weights and returns a look-up-table with the estimated parameters, an estimated parameter table, and cost function values
    
    
    INPUT:
        
        original_formula is a list of formula strings
        formula_weights is a nested list of formula weights
        observable_table is a numpy array with the observables
        observable_table_column_names is a list of strings with names for each of the observables in the observable table
        parameter_position_table is a numpy array matching observable_table in number of rows with the positions of each unknown parameter in its unique table
        parameter_position_table_column_names is a list of strings with the names for each column in parameter_position_table
        individual_parameter_min_max_tables is a list of numpy arrays for each parameter, with min and max values for each parameter
        transfer_function_name is a string with the transfer function to be used
    
    OUTPUT:
        
        estimated_parameters_lut is a minimal-size look-up-table with the estimated parameter values
        estimated_parameter_table is a table matching parameter_position_table in size, where the parameters have been distributed across rows and columns using the estimated_parameters_lut
        cost_function_values is a list of values for each cost function, including the weights
    
    """
    # convert the columnwise indices in "parameter_position_table" to unique indices
    unique_index_table, columnwise_to_unique_index_lut = regularise_indices(parameter_position_table)

    # number of unique parameters to estimate
    length_of_p_vector = columnwise_to_unique_index_lut.shape[0]

    # create a table of min and max parameter values (requires looping through parameter tables and extracting relevant rows)
    # this allows a possible flexible setting of intervals in the future
    p_min_max_table = np.nan * np.zeros((length_of_p_vector, 2))
    for parameter_idx, individual_parameter_min_max_table in enumerate(individual_parameter_min_max_tables):
        current_rows_in_lut = columnwise_to_unique_index_lut[:, 1] == parameter_idx
        current_positions_in_individual_parameter_table = np.int32(
            columnwise_to_unique_index_lut[current_rows_in_lut, 2]
        )
        p_min_max_table[current_rows_in_lut, :] = individual_parameter_min_max_table[
            current_positions_in_individual_parameter_table, :
        ]
    # extract min, max, intiial values
    p_lower = p_min_max_table[:, 0]
    p_upper = p_min_max_table[:, 1]

    # the name of the table that is used in the original_formula string
    # (could be anything as long as the same string is used in swap_names_in_original_formula and cost_function)
    table_name_in_converted_formula = "current_data_table"

    # find rows for which all observable data exist
    # rows_with_all_observables = np.all(~np.isnan(observable_table), axis=1)
    # in this case, the formula should use the first occurence of the same quantity, if it is observed in both "observables" and "parameters"
    converted_formula = swap_names_and_merge_formula(
        original_formula,
        observable_table_column_names,
        parameter_position_table_column_names,
        table_name_in_converted_formula,
        use_observable_if_repeated_and_available=True,
    )

    cost_function_arguments = (
        converted_formula,
        formula_weights,
        observable_table,
        unique_index_table,
        table_name_in_converted_formula,
        lambda x: parameter_transfer_function(x, p_lower, p_upper, False, transfer_function_name),
        True,
    )

    # iterate a few times with different initial values in case some initial value set fails
    max_count = 30
    for counter in range(max_count):

        # creating initial values by randomising
        p_initial = p_lower + np.random.rand(length_of_p_vector) * (p_upper - p_lower)
        # converting to x
        x_initial = parameter_transfer_function(p_initial, p_lower, p_upper, True, transfer_function_name)
        
        # fit the model
        fitted_model = sp.optimize.minimize(cost_function, x_initial, cost_function_arguments, method="BFGS")
        if fitted_model.success:
            p_estimated = parameter_transfer_function(fitted_model.x, p_lower, p_upper, False, transfer_function_name)
            cost_function_values = cost_function(
                parameter_transfer_function(p_estimated, p_lower, p_upper, True, transfer_function_name), *cost_function_arguments[:-1],False
            )
            logging.info(" ... finished with success (cost function values: [%s])." % (', '.join(['%.2f' % (curr_value) for curr_value in cost_function_values])))
            break
        else:
            p_estimated = np.nan * np.zeros(length_of_p_vector)
            cost_function_values = np.nan * np.ones(len(converted_formula))
            if counter < (max_count - 1):
                logging.info(
                    " ... finished with failure (message: {}). Rerunning with different initial values...".format(
                        fitted_model.message
                    )
                )
            else:
                logging.info(
                    " ... finished with failure (message: {}). Entire estimation failed (this should be troubleshot further)...".format(
                        fitted_model.message
                    )
                )
    return (
        np.column_stack((columnwise_to_unique_index_lut, p_initial, p_estimated)),
        p_estimated[unique_index_table],
        cost_function_values,
    )


# %%
def read_and_organise_3d_data(
    current_block_extents,
    block_has_data,
    pixel_axis_north,
    pixel_axis_east,
    stack_info_table,
    stack_info_table_column_names,
    observable_names,
    observable_sources,
    observable_transforms,
    observable_averaging_methods,
    observable_ranges,
    observable_is_required,
    forest_class_sources,
    forest_class_boundaries,
    stack_id_vector,
    forest_class_id_vector,
    space_invariant_parameter_table,
    space_invariant_parameter_names,
    mask_out_area_outside_block=False,
):

    '''
    (
        forest_class_3d,
        observables_3d,
        observable_names,
        space_invariant_parameters_3d,
        space_invariant_parameter_names,
        identifiers_3d,
        identifiers_3d_names,
    ) = read_and_organise_3d_data(
                                current_block_extents,
                                block_has_data,
                                pixel_axis_north,
                                pixel_axis_east,
                                stack_info_table,
                                stack_info_table_column_names,
                                observable_names,
                                observable_sources,
                                observable_transforms,
                                observable_averaging_methods,
                                observable_ranges,
                                observable_is_required,
                                forest_class_sources,
                                forest_class_boundaries,
                                stack_id_vector,
                                forest_class_id_vector,
                                space_invariant_parameter_table,
                                space_invariant_parameter_names,
                                mask_out_area_outside_block=False,
                            )
                                    
    
    
    Uses a newton-based method to fit space variant parameters to observable data and space invariant parameters based on the formula that has to be minimised
        
    INPUT:
        
        current_block_extents contains the east and north boundaries of the current block
        block_has_data contains a bool for each stack indicating whether this stack has any data within the block
        pixel_axis_north, pixel_axis_east are the two axes onto which the read data should be interpolated
        stack_info_table is a table with identifier values for each stack, each stack is represented by a new row
        stack_info_table_column_names is a list of strings with column names for the stack_info_table
        observable_names is a list of strings with names for the observables (quantities to be read)
        observable_sources is a nested list with source data (paths & band_ids) for each stack
        observable_transforms is a list of strings with transforms to be applied to each observable
        observable_averaging_methods is a list of strings with methods of averaging to be applied to each observable
        observable_ranges is a list of accepted ranges for each observable (prior to transformation)
        observable_is_required is a list of bools telling the algorithm if the observable is required for AGB estimation or just optional
        forest_class_sources is a nested list with source data for forest class
        forest_class_boundaries is a list of min and max coordinates for the forest class images
        stack_id_vector is a vector of stack ids matching the space_invariant_parameter_table (below)
        forest_class_id_vector is a vector with forest class ids matching the space_invariant_parameter_table (below)
        space_invariant_parameter_table is a table with the estimated space invariant parameters (which will be rasterised)
        space_invariant_parameter_names is a list of names for the columns of the table above
        mask_out_area_outside_block is a bool indicating whether the area outside block should be set to nan for the rasterised parameter images
        
    OUTPUT:
        
        forest_class_3d is a numpy array with forest class
        observables_3d is a list of numpy_arrays with the observables
        observable_names is a corresponding list of observable names
        space_invariant_parameters_3d is a list of numpy arrays with space invariant parameters (rasterised from the table)
        space_invariant_parameter_names is a corresponding list of parameter names
        identifiers_3d is a numpy array with the identifiers
        identifiers_3d_names is a corresponding list of identifier names
                            
        
    
    
    
    
    
    '''
    
    # create mask for current block
    current_block_mask = np.zeros((len(pixel_axis_north), len(pixel_axis_east)), dtype="bool")
    # set areas within block to true
    current_block_mask[
        np.array(
            [np.where((pixel_axis_north > current_block_extents[3]) & (pixel_axis_north < current_block_extents[2]))[0]]
        ).transpose(),
        np.array(
            [np.where((pixel_axis_east > current_block_extents[0]) & (pixel_axis_east < current_block_extents[1]))[0]]
        ),
    ] = True

    def apply_look_up_table(lut_x, lut_y, output_xs):
        output_y = np.nan * np.zeros(np.prod(output_xs).shape)
        for row_in_lut in range(len(lut_y)):
            positions_in_output = (
                np.prod(
                    [
                        lut_x[row_in_lut, column_in_lut] == output_xs[column_in_lut]
                        for column_in_lut in range(len(output_xs))
                    ]
                )
                == 1
            )
            output_y[positions_in_output] = lut_y[row_in_lut]
        return output_y

    # # derived parameters
    number_of_stacks = stack_info_table.shape[0]
    number_of_observables = len(observable_names)
    # number_of_space_invariant_parameters = len(space_invariant_parameter_names)

    ### READING AND SAMPLING FOREST CLASS DATA
    logging.info("AGB: reading forest class map")

    forest_class_3d = np.zeros((len(pixel_axis_north), len(pixel_axis_east)), dtype="float")

    for file_idx, file_source in enumerate(forest_class_sources):

        forest_class_map_boundaries = forest_class_boundaries[file_idx]

        forest_class_map_is_inside = check_intersection(
            forest_class_map_boundaries[0],
            forest_class_map_boundaries[1],
            forest_class_map_boundaries[2],
            forest_class_map_boundaries[3],
            current_block_extents[0],
            current_block_extents[1],
            current_block_extents[3],
            current_block_extents[2],
        )
        if forest_class_map_is_inside:

            forest_class_3d_curr = np.round(
                interp2d_wrapper(
                    file_source[0], file_source[1] + 1, pixel_axis_east, pixel_axis_north, fill_value=float(0)
                )
            )

            # mean all the fnf tiles
            forest_class_3d = np.ceil(merge_agb_intermediate(forest_class_3d, forest_class_3d_curr, method="nan_mean"))

            logging.info(
                "AGB: reading forest class image data (file {}/{})".format(file_idx + 1, len(forest_class_sources))
            )

            # set all unrealistic values to 0 = non-forest
            forest_class_3d[(forest_class_3d <= 0) | np.isnan(forest_class_3d)] = 0

    forest_class_3d = np.array([forest_class_3d]).transpose([1, 2, 0])

    ### READING OBSERVABLE DATA
    # allocate observable tables (one table in a list for each observable)
    observables_3d = []
    # cycle through observable sources in the stack
    for observable_idx in range(number_of_observables):
        observables_3d.append(np.nan * np.zeros((len(pixel_axis_north), len(pixel_axis_east), number_of_stacks)))

        # number of stacks in current observable (must be either 1 or number_of_stacks)
        current_number_of_stacks = len(observable_sources[observable_idx])

        # check if observable is stacked (i.e., the list length equals number of stacks)
        if current_number_of_stacks == number_of_stacks:

            # cycle over stacks
            for stack_idx in range(number_of_stacks):

                # go ahead only if current stack is (at least partially) contained in the current parameter block:
                if block_has_data[stack_idx]:

                    # cycle over all the equi7 tiles, interpolate over pixel grid and average
                    source_data_interp = np.NaN * np.zeros((len(pixel_axis_north), len(pixel_axis_east)), dtype="float")

                    for file_idx, file_source in enumerate(observable_sources[observable_idx][stack_idx]):

                        source_data_interp_curr = interp2d_wrapper(
                            file_source[0], file_source[1] + 1, pixel_axis_east, pixel_axis_north, fill_value=np.NaN,
                        )

                        source_data_interp = merge_agb_intermediate(
                            source_data_interp, source_data_interp_curr, method="nan_mean"
                        )

                    # masking the stack:
                    source_data_interp[forest_class_3d[:, :, 0] == 0] = np.NaN

                    logging.info(
                        "AGB: reading stacked image data for observable '{}' (stack {}/{}, file {}/{})".format(
                            observable_names[observable_idx],
                            stack_idx + 1,
                            number_of_stacks,
                            file_idx + 1,
                            len(observable_sources[observable_idx][stack_idx]),
                        )
                    )

                    observables_3d[observable_idx][:, :, stack_idx] = transform_function(
                        source_data_interp, observable_ranges[observable_idx], observable_transforms[observable_idx]
                    )

        # otherwise, replicate across stacks
        elif current_number_of_stacks == 1:

            source_data_interp = np.nan * np.zeros((len(pixel_axis_north), len(pixel_axis_east)), dtype="float")

            for file_idx, file_source in enumerate(observable_sources[observable_idx][0]):

                source_data_interp_curr = interp2d_wrapper(
                    file_source[0], file_source[1] + 1, pixel_axis_east, pixel_axis_north, fill_value=np.NaN,
                )

                source_data_interp = merge_agb_intermediate(
                    source_data_interp, source_data_interp_curr, method="nan_mean"
                )

                # masking the stack:
                source_data_interp[forest_class_3d[:, :, 0] == 0] = np.NaN

            logging.info(
                "AGB: reading unstacked image data for observable '{}' (file {}/{})".format(
                    observable_names[observable_idx], file_idx + 1, len(observable_sources[observable_idx][0])
                )
            )
            temporary_transf_image = transform_function(
                source_data_interp, observable_ranges[observable_idx], observable_transforms[observable_idx]
            )
            for stack_idx in range(number_of_stacks):
                observables_3d[observable_idx][:, :, stack_idx] = temporary_transf_image

        # break if this is a required observable and all sampled data are nan
        # (speeds up the reading)
        if np.all(np.isnan(observables_3d[observable_idx])) & observable_is_required[observable_idx]:
            logging.info("AGB: no data for the first required observable, skipping the rest")
            break

    identifiers_3d = []
    for identifier_idx in range(stack_info_table.shape[1]):
        identifiers_3d.append(np.array([[stack_info_table[:, identifier_idx]]]))

    # create maps for the space invariant parameters
    space_invariant_parameters_3d = []
    for parameter_idx, parameter_name in enumerate(space_invariant_parameter_names):
        current_lut = np.column_stack(
            (forest_class_id_vector, stack_id_vector, space_invariant_parameter_table[:, parameter_idx])
        )
        current_lut = np.unique(current_lut, axis=0)
        current_parameter_map_3d = apply_look_up_table(
            current_lut[:, :-1], current_lut[:, -1], (forest_class_3d, identifiers_3d[0])
        )
        if mask_out_area_outside_block:
            current_parameter_map_3d[~current_block_mask] = np.nan
        space_invariant_parameters_3d.append(current_parameter_map_3d)

    # return maps
    return (
        forest_class_3d,
        observables_3d,
        observable_names,
        space_invariant_parameters_3d,
        space_invariant_parameter_names,
        identifiers_3d,
        stack_info_table_column_names,
    )


# %%
def map_space_variant_parameters(
    formula,
    formula_weights,
    forest_class_3d,
    observables_3d,
    observables_3d_names,
    space_invariant_parameters_3d,
    space_invariant_parameters_3d_names,
    identifiers_3d,
    identifiers_3d_names,
    space_variant_parameters_3d_names,
    space_variant_parameters_3d_variabilities,
    space_variant_parameters_3d_limits,
    transfer_function_name,
):
    
    
    '''
    space_variant_parameters_3d, space_variant_parameters_3d_names = map_space_variant_parameters(
                                                                                formula,
                                                                                formula_weights,
                                                                                forest_class_3d,
                                                                                observables_3d,
                                                                                observables_3d_names,
                                                                                space_invariant_parameters_3d,
                                                                                space_invariant_parameters_3d_names,
                                                                                identifiers_3d,
                                                                                identifiers_3d_names,
                                                                                space_variant_parameters_3d_names,
                                                                                space_variant_parameters_3d_variabilities,
                                                                                space_variant_parameters_3d_limits,
                                                                                transfer_function_name,
                                                                            )
                                                                                
    
    
    Uses a newton-based method to fit space variant parameters to observable data and space invariant parameters based on the formula that has to be minimised
        
    INPUT:
        
        formula is a list of strings with the original formulas
        formula_weights is a list or array of floats representing weights for each formula
        forest_class_3d is a numpy array with the forest class data
        observables_3d is a list of numpy arrays with the observable data rasters
        observables_3d_names is a list of strings with names for each observable
        space_invariant_parameters_3d is a list of numpy arrays with the space invariant parameter rasters
        space_invariant_parameters_3d_names is a corresponding list of space invariant parameter names
        identifiers_3d is a list of numpy arrays with identifiers for stacks
        identifiers_3d_names is a corresponding list of identifier names
        space_variant_parameters_3d_names is a list of names for the space variant parameters
        space_variant_parameters_3d_variabilities is a corresponding list of variabilities accross the eight dimesnions
        space_variant_parameters_3d_limits is a corresponding list of 2 element vectors with lower and upper limits for each space variant parameter,
        transfer_function_name is the name of the transfer function to be used for constraining parameters to their limits
        
    OUTPUT:
        
        space_variant_parameters_3d is a list of numpy arrays with the estimated space variant parameter maps
        space_variant_parameters_3d_names is a corresponding list of parameter names
                            
        
    
    
    
    
    
    '''
    
    
    
    def small_change_in_intermediate_parameters_3d(intermediate_parameter, additional_arguments, small_step):
        def cost_function_3d(intermediate_parameter, additional_arguments):
            # print(np.nanmean(intermediate_parameter))
            (
                converted_formula,
                all_observables_3d,
                transfer_function,
                transfer_function_name,
                space_variant_parameter_limits,
                data_list_name,
            ) = additional_arguments
            space_variant_parameter = np.kron(
                np.ones((1, 1, all_observables_3d[0].shape[2])),
                transfer_function(
                    intermediate_parameter, space_variant_parameter_limits[0], space_variant_parameter_limits[1], False,
                    transfer_function_name,
                ),
            )
            data_list = all_observables_3d + [space_variant_parameter]
            return np.sqrt(np.nanmean(eval(converted_formula, {data_list_name: data_list, "np": np}), axis=2))

        cost_function_value = cost_function_3d(intermediate_parameter, additional_arguments)
        cost_function_value_after = cost_function_3d(intermediate_parameter + small_step, additional_arguments)
        cost_function_value_before = cost_function_3d(intermediate_parameter - small_step, additional_arguments)
        small_change = (
            small_step
            * (cost_function_value_after - cost_function_value_before)
            / (2 * (cost_function_value_after - 2 * cost_function_value + cost_function_value_before))
        )
        return np.array([small_change]).transpose([1, 2, 0]), cost_function_value

    if (not len(space_variant_parameters_3d_names) == 1) or (np.any(space_variant_parameters_3d_variabilities[0][1:])):
        logging.error(
            "AGB: map creation is currently not implemented for multiple space-variant parameters or space-variant parameters that change in time or with forest class.",
            exc_info=False,
        )
    else:
        data_list_name = "data_list"
        all_observables_3d_names = observables_3d_names + space_invariant_parameters_3d_names
        converted_formula = swap_names_and_merge_formula_3d(
            formula,
            formula_weights,
            all_observables_3d_names,
            space_variant_parameters_3d_names,
            data_list_name,
            use_observable_if_repeated_and_available=True,
        )
        all_observables_3d = observables_3d + space_invariant_parameters_3d

        # intermediate_parameters_3d = np.pi/4+np.zeros((len(pixel_axis_north),len(pixel_axis_east),1))
        space_variant_parameters_3d_initial = np.mean(space_variant_parameters_3d_limits[0]) * np.ones(
            (observables_3d[0].shape[0], observables_3d[0].shape[1], 1)
        )
        intermediate_parameters_3d = parameter_transfer_function(
            space_variant_parameters_3d_initial,
            space_variant_parameters_3d_limits[0][0],
            space_variant_parameters_3d_limits[0][1],
            True,
            transfer_function_name,
        )

        additional_arguments = (
            converted_formula,
            all_observables_3d,
            parameter_transfer_function,
            transfer_function_name,
            space_variant_parameters_3d_limits[0],
            data_list_name,
        )
        
        
        # this is a tricky part because the small step and max change depend on the quantity that we optimise for
        if transfer_function_name == 'sin2':
                
            small_step = 0.01
            maximal_change_magnitude = 0.03
        else:
            small_step = 0.25
            maximal_change_magnitude = 1
            
            
        number_of_iterations = 1000
        scaling_factor = 0.8
        for ii in np.arange(number_of_iterations):
            small_change, cost_function_value = small_change_in_intermediate_parameters_3d(
                intermediate_parameters_3d, additional_arguments, small_step
            )
            intermediate_parameters_3d = intermediate_parameters_3d - np.maximum(
                -maximal_change_magnitude, np.minimum(maximal_change_magnitude, scaling_factor * small_change)
            )

        space_variant_parameters_3d = parameter_transfer_function(
            intermediate_parameters_3d,
            space_variant_parameters_3d_limits[0][0],
            space_variant_parameters_3d_limits[0][1],
            False,
            transfer_function_name,
        )
        logging.info(
            "AGB: map creation successful (average cost function value: %.2f)" % (np.nanmean(cost_function_value))
        )

        return (
            [space_variant_parameters_3d],
            space_variant_parameters_3d_names,
        )


# %% swap variable names in formulas to elements of a list

  
def swap_names_and_merge_formula_3d(
    original_formulas, weights, observable_names, parameter_names, new_table_name, use_observable_if_repeated_and_available=True
):  
    
    
    '''
    new_formula = swap_names_and_merge_formula_3d(
                    original_formulas, weights, observable_names, parameter_names, new_table_name, use_observable_if_repeated_and_available=True
                )
    
    
    
    Merges individual formula strings and weights to one formula strings and replaces observable and parameter names with elements of a list with pre-defined name.
        
    INPUT:
        
        original_formulas is a list of strings with the original formulas
        weights is a list or array of floats representing weights for each formula
        observable_names is a list of names to formula terms which are to be treated as observables (i.e., known values)
        parameter_names is a list of names to formula terms which are to be treated as parameters (i.e., unknown values)
        new_table_name is a string with the name for the list that contains both the observables and parameters
        use_observable_if_repeated_and_available is a bool indicating whether a name that is observed both in the 
                        observable_names and parameter_names should be treated as observable; by default this is True; 
                        however, this function is going to be depraceted
        
        
    OUTPUT:
        
        new_formula is the output formula as a single string, where each of the original formulas have been squared and summed using the weights
                            
        
    
    
    
    
    
    '''
    original_variable_names = observable_names + parameter_names
    unique_variable_names, name_counts = np.unique(np.array(original_variable_names), return_counts=True)
    new_formula = "0"
    for current_formula, current_weight in zip(original_formulas,weights):
        for unique_variable_name, name_count in zip(unique_variable_names, name_counts):
            if name_count == 1:
                position_in_variable_names_vector = np.where(np.array(original_variable_names) == unique_variable_name)[
                    0
                ][0]
            elif name_count == 2:
                position_in_variable_names_vector = np.where(np.array(original_variable_names) == unique_variable_name)[
                    0
                ][np.int32(~use_observable_if_repeated_and_available)]
            current_formula = current_formula.replace(
                unique_variable_name, new_table_name + ("[%d]" % (position_in_variable_names_vector))
            )
        new_formula = new_formula + " + %f*(%s)**2" % (current_weight,current_formula)
    return new_formula.replace("0 + ", "")


# %%
def match_string_lists(ref_string_list, test_string_list):
    
    
    '''
    match_array = match_string_lists(ref_string_list, test_string_list)
    
    
    
    
    Matches a test string list to a refernece string list element-by-element and returns positions at which the test strings are found in the reference strings
        
    INPUT:
        
        ref_string_list is a list with reference strings
        test_string_list is a list with test strings
        
    OUTPUT:
        
        match_array is a 2D numpy array containing the positions at which strings from test_string_list are found in ref_string_list, -1 is returned if the string is not found
                            
        
    
    
    
    
    
    '''
    
    
    
    is_in = np.zeros((len(ref_string_list), len(test_string_list)))
    for ref_idx, current_ref_string in enumerate(ref_string_list):
        for test_idx, current_test_string in enumerate(test_string_list):
            is_in[ref_idx, test_idx] = current_ref_string.find(current_test_string)
    return is_in



# %% define parameter transfer functions
def parameter_transfer_function(in_vector, p_lower, p_upper, in_vector_is_p=False,transfer_function_name='sin2'):
    
    '''
    out_vector = parameter_transfer_function(in_vector, p_lower, p_upper, in_vector_is_p=False,transfer_function_name='sin2')
    
    
    
    
    Transforms an unconstrained parameter onto a constrained parameter or vice versa.
        
    INPUT:
        
        in_vector is a numpy array with data to be transformed
        p_lower, p_upper are the lower and upper limits for the constrained parameter, they can be either scalars of numpy arrays of the same size as in_vector
        in_vector_is_p determines if the in_vector is the constrained or the unconstrained parameter, if False, then in_vector is unconstrained and out_vector is constrained to [p_lower,p_upper],
                            otherwise, the in_vector is constrained to [p_lower,p_upper] and this function performs the inverted transform
        transfer_function_name is a string defining the transfer function, currently the only allowed transform is "sin2", any other string will return in_vector as out_vector
        
    OUTPUT:
        
        out_vector is the transformed version of in_vector
                            
        
    
    
    
    
    
    '''
     
    
    
    # note: no check of x_vector, p_upper, p_lower is done here,
    # it is assumed that the inputs are correct
    if transfer_function_name == 'sin2':
        if not in_vector_is_p:
            return p_lower + (p_upper - p_lower) * np.sin(in_vector) ** 2
        else:
            return np.arcsin(np.sqrt((in_vector - p_lower) / (p_upper - p_lower)))
    else:
        return in_vector

# %%
def save_human_readable_table(
    path, table, column_names, data_type_lut, table_delimiter, table_precision, table_column_width
):
    
    
    
    '''
    
    save_human_readable_table(
                                path, table, column_names, data_type_lut, table_delimiter, table_precision, table_column_width
                            )
    
    Saves a text file with the table on a nice, human-readable, fixed_width format.
    
    INPUT:
        
        path is a string with the table path
        table is a n x m numpy array with the data
        column_names is a list of length m with the names for each column
        data_type_lut is a list of length 2, where the first contains all possible column names and the second contains associated data types ("f" or "d")
        table_delimiter is the delimiter
        table_precision is the precision for data of type "f"
        table_column_width is the column width
        
    OUTPUT:
        
        no output returned 
        
    
    
    
    '''

    table_format, table_header = get_fmt_and_header(
        column_names, data_type_lut[0], data_type_lut[1], table_delimiter, table_precision, table_column_width
    )

    np.savetxt(path, table, fmt=table_format, delimiter=table_delimiter, header=table_header, comments="")
    path_npy = ".".join(path.split(".")[:-1]) + ".npy"
    np.save(path_npy, table)


# %% function for creating list of format strings and a header for subsequent saving of tables into text files
def get_fmt_and_header(column_names, all_column_groups, all_data_types, delimiter="\t", precision=2, column_width=10):
    
    '''
    
    out_format, out_header = get_fmt_and_header(column_names, all_column_groups, all_data_types, delimiter="\t", precision=2, column_width=10)
    
    Prepares an fmt string and a header string for saving in a nice, human-readable table of fixed column width using np.savetxt
    
    INPUT:
        
        column_names is a list of strings for the current table columns
        all_column_groups is a list of strings for all possible table columns
        all_data_types is a list of strings for each of the columns in all_column_groups (two strings are currently allowed: "d" indicating integer and "f" indicating float)
        delimiter is the delimiter to be used (default is tab)
        precision is the precision to be used for data with type "f" (default is 2)
        column_width is the width of each column (default is 10)
        
        
    OUTPUT:
        
        out_format is the fmt string required by np.savetxt
        out_header is the header required by np.savetxt
        
    
    
    
    '''
    out_format = []
    out_header = []
    for curr_column in column_names:
        for curr_column_group, curr_data_type in zip(all_column_groups, all_data_types):
            if curr_column in curr_column_group:
                if curr_data_type == "d":
                    out_format.append("%s%dd" % (r"%", column_width))
                elif curr_data_type == "f":
                    out_format.append("%s%d.%df" % (r"%", column_width, precision))
                break
        out_header.append("%s%ds" % (r"%", column_width) % curr_column)
    return out_format, delimiter.join(out_header)


# %%


def subset_iterable(iterable_to_subset, validity_mask, return_array=False):
    '''
    
    out_iterable = subset_iterable(in_iterable,validity_mask, return_array=False)
    
    Applies validity mask to an interable, where the iterable can be numpy vector, a list, or a tuple
    
    INPUT:
        
        in_iterable is any iterable of length N (numpy vector, list, tuple)
        validity_mask is a boolean vector of length with True marking that this particular item has to be taken out
        return_array is a boolean indicating if the returned iterable should be a numpy vector
        
        
    OUTPUT:
        
        out_iterable is an iterable of length=sum(validity_mask) containing a subset of the in_iterable matching the validity_mask, if return_array=True then out_iterable is a numpy array
        
    
    
    
    '''
    out = [value for value, flag in zip(iterable_to_subset, validity_mask) if flag]
    if return_array:
        return np.array(out)
    else:
        return out


# %%


## in the future, improve this so it can handle polygons etc
def check_block_for_data_and_cal(block_extents, stack_boundaries, calibration_boundaries):
    '''
    
    block_has_data, block_has_cal = check_block_for_data_and_cal(block_extents, stack_boundaries, calibration_boundaries)
    
    Checks if the current block has any stack data and calibration data
    
    INPUT:
        
        block_extents is a list with min and max easting and northing coordinates
        stack_boundaries is an array with min and max easting and northing coordinates for each stack
        calibration_boundaries is an array with min and max easting and northing coordinates for each calibration dataset
        
    OUTPUT:
        
        block_has_data is a vector of bools with length equal to the number of stacks, indicating if the block contains each stack
        block_has_cal is a vector of bools with length equal to the number of calibration datasets, indicating if the block contains the calibration dataset
    
    
    
    
    '''

    # cycle through stacks and check that there are some data within current block
    block_has_data = np.zeros(stack_boundaries.shape[0], dtype="bool")
    for stack_idx, stack_boundary in enumerate(stack_boundaries):

        # go ahead only if current parameter block is at least partially contained in the data stack:
        block_has_data[stack_idx] = check_intersection(
            stack_boundary[0],
            stack_boundary[1],
            stack_boundary[2],
            stack_boundary[3],
            block_extents[0],
            block_extents[1],
            block_extents[3],
            block_extents[2],
        )

    # cycle through cals and see that there are some cals within current block
    block_has_cal = np.zeros(calibration_boundaries.shape[0], dtype="bool")
    for calibration_idx, calibration_boundary in enumerate(calibration_boundaries):

        # go ahead only if current parameter block is at least partially contained in the data stack:
        block_has_cal[calibration_idx] = check_intersection(
            calibration_boundary[0],
            calibration_boundary[1],
            calibration_boundary[2],
            calibration_boundary[3],
            block_extents[0],
            block_extents[1],
            block_extents[3],
            block_extents[2],
        )

    return block_has_data, block_has_cal


# %%
def compute_block_processing_order(
    block_corner_coordinates_east,
    block_corner_coordinates_north,
    block_size_east,
    block_size_north,
    calibration_area_coordinates,
    stack_data_coordinates,
):
    
    '''
    
    block_order = compute_block_processing_order(
                                                    block_corner_coordinates_east,
                                                    block_corner_coordinates_north,
                                                    block_size_east,
                                                    block_size_north,
                                                    calibration_area_coordinates,
                                                    stack_data_coordinates,
                                                )
    
    Calculates the order in which blocks are to be processed, based on the available calibration area coordinates and stack data coordinates
    
    INPUT:
        
        block_corner_coordinates_east,block_corner_coordinates_north are vectors of matching size specifying the corner coordinates of each block (in metres)
        block_size_east,block_size_north are scalars specifying the block size in each direction (in metres)
        calibration_area_coordinates, stack_data_coordiantes are arrays with coordinates for calibration and stack data, respectively, with one row for each calibration dataset/stack
        
    OUTPUT:
        
        block_order is a vector of the same length as the number of blocks containing block IDs in the order in which they should be run
    
    
    
    
    '''
    # in the future, this should be capable of reading polygons for both calibration areas and stack data
    # now, it's just a wrapper around an old function
    (current_block_index, block_order) = compute_processing_blocs_order(
        calibration_area_coordinates,
        block_corner_coordinates_east,
        block_size_east,
        block_corner_coordinates_north,
        block_size_north,
    )
    return block_order


# %% transform functions
def transform_function(in_data, interval, kind, do_forward=True):
    
    
    '''
    
    out_data = transform_function(in_data,interval,kind,do_forward=True)
    
    
    Transforms data within a certain interval using a transformation
    
    
    INPUT:
    in_data is a scalar or a numpy array with the data to be transformed
    interval is a list/vector/tuple of length 2 with the lower and upper interval bounds (prior to transformation); outside this interval, the transformed values are nan
    kind is a string with the type of transformation
    do_forward is a bool determining the direction of the transformation (True means forward transformation, False means inverse transformation)
    
    OUTPUT:
    out_data is the transformed dataset with the same shape as in_data and with nans for values outside the interval
    
    
    
    
    
    
    
    '''
    out_data = np.copy(in_data)
    out_data[(out_data < interval[0]) | (out_data > interval[1])] = np.nan
    # note: no check of data and kind is done here,
    # it is assumed that the inputs are correct
    if kind == "db":
        if do_forward:
            out_data = 10 * np.log10(out_data)
        else:
            out_data = 10 ** (0.1 * in_data)
    elif kind == "-db":
        if do_forward:
            out_data = -10 * np.log10(in_data)
        else:
            out_data = 10 ** (-0.1 * in_data)
    elif kind == "-2db":
        if do_forward:
            out_data = -10 * np.log10(2 * in_data)
        else:
            out_data = 10 ** (-0.1 * in_data) / 2
    elif kind == "cosdb":
        if do_forward:
            out_data = 10 * np.log10(np.cos(in_data))
        else:
            out_data = np.arccos(10 ** (0.1 * in_data))
    elif kind == "cos":
        if do_forward:
            out_data = np.cos(in_data)
        else:
            out_data = np.arccos(in_data)
    elif kind == "round":
        out_data = np.round(in_data)
        
    else:
        out_data = np.copy(in_data)
    return out_data


# %% function for calculating a given statistic for polygons of arbitrary shape
def stats_on_polygons(data, pixel_axis_east, pixel_axis_north, reference_polygons, method):
    
    '''
    
    stats = stats_on_polygons(data_image,pixel_axis_east,pixel_axis_north,polygons,method)
    
    
    Calculates a statistic from an image over all polygons 
    
    
    INPUT:
    data_image is a 2D array
    pixel_axis_east, pixel_axis_north are two 1D arrays defining the east and north axes for the data_image
    polygons is a list of N Polygon objects over which the statistic is calculated
    method is a string with the statistic to be calculated (see "stats_on_polygons" to see the statistics that can be calculated)
    
    OUTPUT:
    stats isa 1D array with length N containing the statistic calculated for each polygon
    
    
    
    
    
    
    
    '''

    # initialize inputs:
    Nx, Ny = data.shape

    data_east_min = min(pixel_axis_east.flatten())
    data_east_delta = (max(pixel_axis_east.flatten()) - min(pixel_axis_east.flatten())) / Nx
    data_north_in = pixel_axis_north.flatten()[0]
    data_north_delta = (pixel_axis_north.flatten()[-1] - pixel_axis_north.flatten()[0]) / Ny
    # input data geotransform:
    data_geotransform = [data_east_min, data_east_delta, 0, data_north_in, 0, data_north_delta]
    # Setup working spatial reference
    sr_wkt = 'LOCAL_CS["arbitrary"]'
    sr = osr.SpatialReference(sr_wkt)

    # initiale output stats vector
    polygon_means_vec = np.nan * np.zeros(len(reference_polygons))

    for index, polygon in enumerate(reference_polygons):
        
        # Create a memory raster (gdal raster) to rasterize into.
        raster_support_driver = gdal.GetDriverByName("MEM").Create("", Ny, Nx, 1, gdal.GDT_Byte)
        raster_support_driver.SetGeoTransform(data_geotransform)
        raster_support_driver.SetProjection(sr_wkt)

        # Create a memory ploygon layer (ogr vector) to rasterize from.
        poly_layer_support_driver = ogr.GetDriverByName("Memory").CreateDataSource("wrk")
        poly_layer = poly_layer_support_driver.CreateLayer("poly", srs=sr)

        # Add a polygon to the layer.
        if isinstance(polygon, str):
            wkt_geom = polygon
        elif isinstance(polygon, Polygon):
            wkt_geom = polygon.wkt
        feat = ogr.Feature(poly_layer.GetLayerDefn())
        feat.SetGeometryDirectly(ogr.Geometry(wkt=wkt_geom))
        poly_layer.CreateFeature(feat)

        # Run the rasterization of polygon over raster data driver
        gdal.RasterizeLayer(raster_support_driver, [1], poly_layer, burn_values=[1])
        # read the mask from the rasterized polygon over data driver
        bandmask = raster_support_driver.GetRasterBand(1)
        datamask = bandmask.ReadAsArray().astype(np.bool)
        if np.any(datamask):
            # Calculate statistics of zonal raster
            if method == "mean":
                polygon_means_vec[index] = np.mean(data[datamask])
    
            elif method == "nan_mean":
                polygon_means_vec[index] = np.nanmean(data[datamask])
    
            elif method == "median":
                polygon_means_vec[index] = np.median(data[datamask])
    
            elif method == "nan_median":
                polygon_means_vec[index] = np.nanmedian(data[datamask])
            elif method == "mode":
                polygon_means_vec[index] = sp.stats.mode(data[datamask])[0]
            elif method == "mode_knp":
                if np.any(data[datamask]<=0):
                    polygon_means_vec[index] = np.nan
                else:
                    polygon_means_vec[index] = sp.stats.mode(data[datamask])[0]
                    
            elif method == "sum":
                polygon_means_vec[index] = np.sum(data[datamask])
            elif method == "range":
                polygon_means_vec[index] = np.max(data[datamask])-np.min(data[datamask])

        raster_support_driver = None
        poly_layer_support_driver = None

    return polygon_means_vec


# %% function for calculating a given statistic on all samples on a grid and all polygons
def stats_on_all_samples(
    data_image,
    data_resolution,
    pixel_axis_east,
    pixel_axis_north,
    polygons,
    method,
):
    '''
    
    stats = stats_on_all_samples(data_image,data_resolution,pixel_axis_east,pixel_axis_north,polygons,method)
    
    
    Calculates a statistic from an image over all polygons with size larger or equal to the resolution
    
    
    INPUT:
    data_image is a 2D array
    data_resolution is resolution in metres
    pixel_axis_east, pixel_axis_north are two 1D arrays defining the east and north axes for the data_image
    polygons is a list of N Polygon objects over which the statistic is calculated
    method is a string with the statistic to be calculated (see "stats_on_polygons" to see the statistics that can be calculated)
    
    OUTPUT:
    stats isa 1D array with length N containing the statistic calculated for each polygon (nan is returned if polygon.area < data_resolution**2)
    
    
    
    
    
    
    
    '''
    stats_polygons = stats_on_polygons(data_image, pixel_axis_east, pixel_axis_north, polygons, method)
    
    def get_polygon_areas(polygons):
        return np.array([polygon.area for polygon in polygons])
    # kill the polygons with area smaller than the resolution cell 
    stats_polygons[get_polygon_areas(polygons)<data_resolution**2] = np.nan        
    
    return stats_polygons

