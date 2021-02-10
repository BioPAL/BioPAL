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
    mean_on_rois,
    check_intersection,
    interp2d_wrapper,
    merge_agb_intermediate,
    compute_processing_blocs_order,
)


# %% reading tables
def sample_and_tabulate_data(
    block_extents,  # extent of the current area for which the table is created
    pixel_axis_east,  # east north axes onto which data are interpolated
    pixel_axis_north,
    sample_axis_east,  # east north axes of the regular sampling grid
    sample_axis_north,
    sample_size_east,  # east north extents of the samples
    sample_size_north,
    additional_sampling_polygons,  # additional arbitrarily shaped polygons
    stack_in_block,  # flags whether each stack is in current block
    stack_info_table,  # info table with stack properties (stack id, headings, etc., defining the acquisition parameters)
    stack_info_table_column_names,  # column names for the abovementioned table
    forest_class_boundaries,
    forest_class_sources,
    observable_names,  # observable names in formula
    observable_is_required,  # flags whether the observable must exist in each table row
    observable_sources,  # paths to equi7 tif files paired with band IDs
    observable_transforms,  # transform function to apply to observable
    observable_averaging_methods,  # averaging method (most commonly 'mean', but could be other if required (e.g., for slope aspect angle))
    observable_ranges,  # permitted ranges, outside those the observable is set to nan
    parameter_names,  # parameter names in formula
    parameter_limits,  # permissible parameter intervals
    parameter_variabilities,  # parameter variabilities across all dimensions
    number_of_subsets,  # number of subsets to use (used to allocate columns in parameter tables)
):

    # derived parameters
    number_of_stacks = stack_info_table.shape[0]
    stack_id_vector = stack_info_table[:, 0]
    number_of_samples = len(sample_axis_east) * len(sample_axis_north) + len(additional_sampling_polygons)
    number_of_parameters = len(parameter_names)
    number_of_observables = len(observable_names)
    pixel_mesh_east, pixel_mesh_north = np.meshgrid(pixel_axis_east, pixel_axis_north)

    # defining names for identifiers (sampleID & forest class ID and then all columns from stack info table)
    identifier_names = ["sample_id", "forest_class_id"] + stack_info_table_column_names
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

    ### READING AND SAMPLING FOREST CLASS DATA
    logging.info("AGB: sampling forest class data")

    forest_class_map_interp = np.zeros((len(pixel_axis_north), len(pixel_axis_east)), dtype="float")

    for file_idx, file_source in enumerate(forest_class_sources):

        forest_class_map_boundaries = forest_class_boundaries[file_idx]

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
                    file_source[0], file_source[1] + 1, pixel_axis_east, pixel_axis_north, fill_value=float(0)
                )
            )

            # mean all the fnf tiles
            forest_class_map_interp = np.ceil(
                merge_agb_intermediate(forest_class_map_interp, forest_class_map_interp_curr, method="nan_mean")
            )

            logging.info("AGB: sampling forest class data (file {}/{})".format(file_idx + 1, len(forest_class_sources)))

            # set all unrealistic values to 0 = non-forest
            forest_class_map_interp[(forest_class_map_interp <= 0) | np.isnan(forest_class_map_interp)] = 0

        # err_str = 'Cannot find any forest class map falling in current block coordinates.'
        # logging.error(err_str)
        # raise ImportError(err_str)

    # sampling forest class map
    temp_forest_class_vector = np.int32(
        np.round(
            stats_on_all_samples(
                forest_class_map_interp,
                pixel_mesh_east,
                pixel_mesh_north,
                sample_axis_east,
                sample_size_east,
                sample_axis_north,
                sample_size_north,
                additional_sampling_polygons,
                "mode",
            )
        )
    )
    # repeating it across stacks and inserting into observable data table (note: forest class assumed constant across stacks)
    identifier_table[:, 1] = np.kron(temp_forest_class_vector, np.ones(number_of_stacks))

    ### READING OBSERVABLE DATA
    # allocate observable table
    observable_table = np.nan * np.zeros((number_of_samples * number_of_stacks, number_of_observables))
    # cycle through observable sources in the stack
    for observable_idx in range(number_of_observables):
        # number of stacks in current observable (must be either 1 or number_of_stacks)
        current_number_of_stacks = len(observable_sources[observable_idx])

        # check if observable is stacked (i.e., the list length equals number of stacks)
        if current_number_of_stacks == number_of_stacks:

            # cycle over stacks
            for stack_idx in range(number_of_stacks):

                # go ahead only if current stack is (at least partially) contained in the current parameter block:
                if stack_in_block[stack_idx]:

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
                    source_data_interp[forest_class_map_interp == 0] = np.NaN

                    logging.info(
                        "AGB: sampling and transforming stacked data for observable '{}' (stack {}/{}, file {}/{})".format(
                            observable_names[observable_idx],
                            stack_idx + 1,
                            number_of_stacks,
                            file_idx + 1,
                            len(observable_sources[observable_idx][stack_idx]),
                        )
                    )

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
                        observable_transforms[observable_idx],
                    )
                    # find rows where the stack id matches the current stack id
                    current_rows = identifier_table[:, 2] == stack_id_vector[stack_idx]
                    # fill out the table
                    observable_table[current_rows, observable_idx] = temp_transformed_sampled_data

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

                logging.info(
                    "AGB: sampling and transforming unstacked data for observable '{}' (file {}/{})".format(
                        observable_names[observable_idx], file_idx + 1, len(observable_sources[observable_idx][0])
                    )
                )

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
                    method="nan_mean",
                ),
                observable_ranges[observable_idx],
                observable_transforms[observable_idx],
            )
            # fill out the table
            observable_table[:, observable_idx] = np.kron(temp_transformed_sampled_data, np.ones(number_of_stacks))

        # break if this is a required observable and all sampled data are nan
        # (speeds up the reading)
        if np.all(np.isnan(observable_table[:, observable_idx])) & observable_is_required[observable_idx]:
            logging.info("AGB: no data for the first required observable, skipping the rest")
            break

    # observable_table = np.nan*observable_table # this is just a dummy thing to check the behaviour of this function in case of lack of data

    # mark rows in observable data table that have negative identifiers, nan-valued sar observables, infinite sar observables, or negative agb values
    invalid_rows = (
        np.any(identifier_table < 0, axis=1)
        | np.any(np.isnan(observable_table[:, observable_is_required]), axis=1)
        | np.any(~np.isfinite(observable_table[:, observable_is_required]), axis=1)
    )
    # exclude invalid rows
    observable_table = observable_table[~invalid_rows, :]
    identifier_table = identifier_table[~invalid_rows, :]
    # number of rows in data table
    number_of_rows_in_observable_table = observable_table.shape[0]

    ### PREPARING PARAMETER TABLES
    parameter_property_names = ["lower_limit", "upper_limit", "initial_value"] + [
        "estimate_%d" % (ii) for ii in np.arange(number_of_subsets) + 1
    ]
    parameter_position_names = ["row_" + parameter_name for parameter_name in parameter_names]
    parameter_tables = []
    parameter_table_columns = []
    parameter_position_table = np.nan * np.zeros((number_of_rows_in_observable_table, number_of_parameters))
    # creating parameter matrices
    for parameter_idx, parameter_variability in enumerate(parameter_variabilities):
        # take out only the relevant identifiers (the ones that change as per parameter variability)
        # and create column names by adding four additional columns: min, max and initial value, and estimated value (later set to NaN)
        parameter_table_columns.append(
            np.concatenate((np.array(identifier_names)[parameter_variability], np.array(parameter_property_names)))
        )
        # create the minimal ID table (without unnecessary columns for those dimension across which the parameter doesn't change)
        temp_ids_table = identifier_table[:, np.where(parameter_variability)[0]]
        # create the last four columns
        temp_minmax_table = np.array([parameter_limits[parameter_idx]]) * np.ones(
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
        observable_names,
        identifier_table,
        identifier_names,
        parameter_position_table,
        parameter_position_names,
        parameter_tables,
        parameter_table_columns,
    )


# %%


def fit_formula_to_random_subsets(
    formula,
    formula_weights,
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
):

    ### CREATE CALIBRATION AND ESTIMATION SUBSETS
    logging.info("AGB: creating {} subsets".format(number_of_subsets))

    # select rows with available agb information as calibration data and those without as estimation data
    calibration_rows = np.where(np.all(~np.isnan(observable_table), axis=1))[0]
    estimation_rows = np.where(np.any(np.isnan(observable_table), axis=1))[0]
    calibration_sample_ids = np.unique(identifier_table[calibration_rows, 0])
    estimation_sample_ids = np.unique(identifier_table[estimation_rows, 0])

    # calculate subset sizes
    # estimation_subset_size = np.int32(np.ceil(len(estimation_sample_ids)/100*proc_conf.AGB.fraction_of_roi_per_test))
    # calibration_subset_size = np.int32(np.ceil(len(calibration_sample_ids)/100*proc_conf.AGB.fraction_of_cal_per_test))
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
            logging.info(
                "AGB: warning: number of accepted subsets ({}) is less than the required number of subsets ({}).".format(
                    number_of_accepted_subsets, number_of_subsets
                )
            )

        ### FIRST, WE PERFORM THE ESTIMATION OF PARAMETERS FOR SUBSETS

        # find observables and parameters in formula and create
        # vectors for selecting parameters and observables that exist in formula
        observables_in_formula = np.any(match_string_lists(formula, observable_names) >= 0, axis=0)
        observables_in_parameters = np.any(match_string_lists(parameter_names, observable_names) >= 0, axis=0)
        parameters_in_formula = np.any(match_string_lists(formula, parameter_names) >= 0, axis=0)
        # find parameters that do not change between samples
        space_invariant_parameters = False == np.column_stack(parameter_variabilities)[0, :]
        # some more flag vectors
        space_invariant_parameters_in_formula = space_invariant_parameters & parameters_in_formula
        space_variant_parameters_in_formula = ~space_invariant_parameters & parameters_in_formula
        observables_in_formula_not_in_parameters = observables_in_formula & ~observables_in_parameters

        # loop through calibration subsets
        for subset_idx, current_subset in enumerate(subset_indexing_vectors):

            logging.info(
                "AGB: minimising cost function for subset {} out of {}...".format(subset_idx + 1, number_of_subsets)
            )

            # subset parameter and observable tables
            current_parameter_position_table = parameter_position_table[current_subset, :][:, parameters_in_formula]
            current_parameter_names = subset_iterable(parameter_names, parameters_in_formula, False)
            current_observable_table = observable_table[current_subset, :][:, observables_in_formula]
            current_observable_names = subset_iterable(observable_names, observables_in_formula, False)
            # create min-max tables
            individual_parameter_min_max_tables = subset_iterable(parameter_tables, parameters_in_formula, False)
            for parameter_idx in range(len(individual_parameter_min_max_tables)):
                individual_parameter_min_max_tables[parameter_idx] = individual_parameter_min_max_tables[parameter_idx][
                    :, -number_of_subsets - 3 : -number_of_subsets - 1
                ]

            # estimate both parameters and AGB for the subset
            (current_lut_all_parameters, _, _,) = fit_formula_to_table_data(
                formula,
                formula_weights,
                current_observable_table,
                current_observable_names,
                current_parameter_position_table,
                current_parameter_names,
                individual_parameter_min_max_tables,
            )

            # fill out parameter tables with estimates of space invariant parameters
            for current_parameter_idx in np.where(space_invariant_parameters_in_formula)[0]:
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

        # loop through calibration subsets
        for subset_idx in range(number_of_subsets):

            logging.info(
                "AGB: estimating space-variant parameters for all samples using space-invariant parameter set {} out of {}...".format(
                    subset_idx + 1, number_of_subsets
                )
            )

            # create a table with all space invariant parameters
            space_invariant_parameter_names = subset_iterable(
                parameter_names, space_invariant_parameters_in_formula, False
            )
            current_space_invariant_parameter_table = []
            for current_parameter_idx in np.where(space_invariant_parameters_in_formula)[0]:
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
                (observable_table[:, observables_in_formula_not_in_parameters], current_space_invariant_parameter_table)
            )
            new_observable_names = (
                subset_iterable(observable_names, observables_in_formula_not_in_parameters, False)
                + space_invariant_parameter_names
            )

            # now, only the space-variant parameters are treated us unknown parameters
            # the corresponding tables and lists are now created
            new_parameter_position_table = parameter_position_table[:, space_variant_parameters_in_formula]
            new_parameter_names = subset_iterable(parameter_names, space_variant_parameters_in_formula, False)
            new_individual_parameter_min_max_tables = subset_iterable(
                parameter_tables, space_variant_parameters_in_formula, False
            )
            for parameter_idx in range(len(new_individual_parameter_min_max_tables)):
                new_individual_parameter_min_max_tables[parameter_idx] = new_individual_parameter_min_max_tables[
                    parameter_idx
                ][:, -number_of_subsets - 3 : -number_of_subsets - 1]

            # estimate space variant parameters for all samples
            (current_lut_space_variant_parameters, space_variant_parameter_table, _,) = fit_formula_to_table_data(
                formula,
                formula_weights,
                new_observable_table,
                new_observable_names,
                new_parameter_position_table,
                new_parameter_names,
                new_individual_parameter_min_max_tables,
            )

            # fill out parameter tables with estimates of space invariant parameters
            for current_parameter_idx in np.where(space_variant_parameters_in_formula)[0]:
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

        logging.info("AGB: estimating space-invariant parameters using space-variant parameter estimate")

        # names for table columns
        space_variant_parameter_names = subset_iterable(parameter_names, space_variant_parameters_in_formula, False)
        # create table with all space invariant parameters
        current_space_variant_parameter_table = []
        for position_in_parameter_table_list in np.where(space_variant_parameters_in_formula)[0]:
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
            (observable_table[:, observables_in_formula_not_in_parameters], current_space_variant_parameter_table)
        )
        new_observable_names = (
            subset_iterable(observable_names, observables_in_formula_not_in_parameters, False)
            + space_variant_parameter_names
        )

        new_parameter_position_table = parameter_position_table[:, space_invariant_parameters_in_formula]
        new_parameter_names = subset_iterable(parameter_names, space_invariant_parameters_in_formula, False)
        new_individual_parameter_min_max_tables = subset_iterable(
            parameter_tables, space_invariant_parameters_in_formula, False
        )
        for parameter_idx in range(len(new_individual_parameter_min_max_tables)):
            new_individual_parameter_min_max_tables[parameter_idx] = new_individual_parameter_min_max_tables[
                parameter_idx
            ][:, -number_of_subsets - 3 : -number_of_subsets - 1]

        # estimate space variant parameters for all samples
        (_, space_invariant_parameter_table, _,) = fit_formula_to_table_data(
            formula,
            formula_weights,
            new_observable_table,
            new_observable_names,
            new_parameter_position_table,
            new_parameter_names,
            new_individual_parameter_min_max_tables,
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


# function for converting columnwise indices (which are repeated within the same column if they represent identical values,
# but which may be repeated across different columns without meaning that they represent identical values)
# to unique indices (which are only repeated within the same table if they are to have identical values)
def regularise_indices(columnwise_index_table):

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


def cost_function(
    x_vector, converted_formulas, formula_weights, observable_table, index_table, name_of_table_in_converted_formula, transfer_function, return_one_value=True
):
    
    
    if return_one_value:
        
        p_vector = transfer_function(x_vector)
        table_in_converted_formula = np.column_stack((observable_table, p_vector[index_table]))
        final_expression = '0'
        for converted_formula,formula_weight in zip(converted_formulas,formula_weights):
            final_expression += '+%.18f*np.nanmean((%s)**2)' % (formula_weight/np.sum(formula_weights),converted_formula)
        final_expression = 'np.sqrt((%s))' % (final_expression)
        total_cost = eval(final_expression, {name_of_table_in_converted_formula: table_in_converted_formula, "np": np})
        return total_cost
    
    else:
        
        p_vector = transfer_function(x_vector)
        table_in_converted_formula = np.column_stack((observable_table, p_vector[index_table]))
        final_expression = '0'
        for converted_formula,formula_weight in zip(converted_formulas,formula_weights):
            final_expression += ',%.18f*np.nanmean((%s)**2)' % (formula_weight/np.sum(formula_weights)*len(formula_weights),converted_formula)
        final_expression = 'np.sqrt(np.array([%s]))' % (final_expression)
        total_cost = eval(final_expression, {name_of_table_in_converted_formula: table_in_converted_formula, "np": np})[1:]
        return total_cost


def fit_formula_to_table_data(
    original_formula,
    formula_weights,
    observable_table,
    observable_table_column_names,
    parameter_position_table,
    parameter_position_table_column_names,
    individual_parameter_min_max_tables,
):
   

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

    # if np.all(rows_with_all_observables):
    # if all rows have all data, the total cost function uses only the calibration formula
    cost_function_arguments = (
        converted_formula,
        formula_weights,
        observable_table,
        unique_index_table,
        table_name_in_converted_formula,
        lambda x: parameter_transfer_function(x, p_lower, p_upper, False),
        True,
    )
    # else:
    #     # add an estimation ssd
    #     converted_estimation_formula = swap_names_and_merge_formula(
    #         original_formula,
    #         observable_table_column_names,
    #         parameter_position_table_column_names,
    #         table_name_in_converted_formula,
    #         use_observable_if_repeated_and_available=False,
    #     )
    #     cost_function_arguments = (
    #         [converted_calibration_formula, converted_estimation_formula],
    #         [observable_table[rows_with_all_observables, :], observable_table],
    #         [unique_index_table[rows_with_all_observables, :], unique_index_table],
    #         table_name_in_converted_formula,
    #         lambda x: parameter_transfer_function(x, p_lower, p_upper, False),
    #     )

    # iterate a few times with different initial values in case some initial value set fails
    max_count = 10
    for counter in range(max_count):

        # creating initial values by randomising
        p_initial = p_lower + np.random.rand(length_of_p_vector) * (p_upper - p_lower)
        # converting to x
        x_initial = parameter_transfer_function(p_initial, p_lower, p_upper, True)
        
        # fit the model
        fitted_model = sp.optimize.minimize(cost_function, x_initial, cost_function_arguments, method="BFGS")
        if fitted_model.success:
            p_estimated = parameter_transfer_function(fitted_model.x, p_lower, p_upper, False)
            cost_function_values = cost_function(
                parameter_transfer_function(p_estimated, p_lower, p_upper, True), *cost_function_arguments[:-1],False
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
):
    def small_change_in_intermediate_parameters_3d(intermediate_parameter, additional_arguments, small_step):
        def cost_function_3d(intermediate_parameter, additional_arguments):
            # print(np.nanmean(intermediate_parameter))
            (
                converted_formula,
                all_observables_3d,
                transfer_function,
                space_variant_parameter_limits,
                data_list_name,
            ) = additional_arguments
            space_variant_parameter = np.kron(
                np.ones((1, 1, all_observables_3d[0].shape[2])),
                transfer_function(
                    intermediate_parameter, space_variant_parameter_limits[0], space_variant_parameter_limits[1], False
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
        )

        additional_arguments = (
            converted_formula,
            all_observables_3d,
            parameter_transfer_function,
            space_variant_parameters_3d_limits[0],
            data_list_name,
        )
        small_step = 0.001
        maximal_change_magnitude = 0.03
        number_of_iterations = 100
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
        )
        logging.info(
            "AGB: map creation successful (average cost function value: %.2f)" % (np.nanmean(cost_function_value))
        )

        return (
            [space_variant_parameters_3d],
            space_variant_parameters_3d_names,
        )


#% swap variable names in formulas to slices of an array
def swap_names_and_merge_formula_3d(
    original_formulas, observable_names, parameter_names, new_table_name, use_observable_if_repeated_and_available=True
):
    original_variable_names = observable_names + parameter_names
    unique_variable_names, name_counts = np.unique(np.array(original_variable_names), return_counts=True)
    new_formula = "0"
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
                unique_variable_name, new_table_name + ("[%d]" % (position_in_variable_names_vector))
            )
        new_formula = new_formula + " + (%s)**2" % (current_formula)
    return new_formula.replace("0 + ", "")


# %%
def match_string_lists(ref_string_list, test_string_list):
    is_in = np.zeros((len(ref_string_list), len(test_string_list)))
    for ref_idx, current_ref_string in enumerate(ref_string_list):
        for test_idx, current_test_string in enumerate(test_string_list):
            is_in[ref_idx, test_idx] = current_ref_string.find(current_test_string)
    return is_in


# # %% cost function from formulas and tables
# def cost_function(x_vector,converted_formulas,observable_tables,index_tables,name_of_table_in_converted_formula,transfer_function):
#     # number of independent formulas
#     number_of_formulas = len(converted_formulas)
#     # convert unconstrained x vector to constrained p vector
#     p_vector = transfer_function(x_vector)
#     # allocate total cost
#     total_cost = 0
#     # loop through different formulas and associated tables
#     for observable_table,index_table,converted_formula in zip(observable_tables,index_tables,converted_formulas):
#         table_in_converted_formula = np.column_stack((observable_table,p_vector[index_table]))
#         # evaluate the cost function formula, average across different measurements and add to the total cost
#         total_cost += np.mean(eval(converted_formula,{name_of_table_in_converted_formula:table_in_converted_formula,'np':np}))
#         # print(eval(converted_formula,{name_of_table_in_converted_formula:table_in_converted_formula}))
#     # return weighed with the number of formulas
#     return np.sqrt(total_cost/number_of_formulas)
# #%% function for converting columnwise indices (which are repeated within the same column if they represent identical values,
# # but which may be repeated across different columns without meaning that they represent identical values)
# # to unique indices (which are only repeated within the same table if they are to have identical values)
# def regularise_indices(columnwise_index_table):

#     offset = 0
#     unique_index_table = [] # position in a single x-vector
#     columnwise_to_unique_index_lut = [] # lut for converting between parameter id and column and parameter position in x
#     for column_idx,parameter_column in enumerate(columnwise_index_table.transpose()):
#         # add -1 at the beginning to avoid vectors with one element (which will not work with interp1d)
#         old_indices = np.concatenate((-1*np.ones(1),np.unique(parameter_column)))
#         # new indices is a simple sequence from 0 to number of parameters-1 + offset due to previous parameters
#         new_indices = np.arange(len(old_indices))-1+offset
#         # convert parameter indices and add to the list
#         unique_index_table.append(sp.interpolate.interp1d(old_indices,new_indices,kind='nearest')(parameter_column))
#         # save the lut, removing the first, unnecessary element
#         columnwise_to_unique_index_lut.append(
#             np.column_stack((new_indices[old_indices>-1],
#                              column_idx+np.zeros(len(old_indices[old_indices>-1])),
#                              old_indices[old_indices>-1],
#                              )))
#         # update offset based on current parameter column
#         offset = np.max(new_indices)+1
#     # convert the list of vectors to an array
#     unique_index_table = np.int32(np.column_stack(unique_index_table))
#     # stack all luts to one lut
#     columnwise_to_unique_index_lut = np.row_stack(columnwise_to_unique_index_lut)
#     # # length of beta vector
#     # length_of_p_vector = columnwise_to_unique_index_lut.shape[0]
#     return unique_index_table,columnwise_to_unique_index_lut


# %% define parameter transfer functions
def parameter_transfer_function(in_vector, p_lower, p_upper, in_vector_is_p=False):
    # note: no check of x_vector, p_upper, p_lower is done here,
    # it is assumed that the inputs are correct
    if not in_vector_is_p:
        return p_lower + (p_upper - p_lower) * np.sin(in_vector) ** 2
    else:
        return np.arcsin(np.sqrt((in_vector - p_lower) / (p_upper - p_lower)))


# %%
def save_human_readable_table(
    path, table, column_names, data_type_lut, table_delimiter, table_precision, table_column_width
):

    table_format, table_header = get_fmt_and_header(
        column_names, data_type_lut[0], data_type_lut[1], table_delimiter, table_precision, table_column_width
    )

    np.savetxt(path, table, fmt=table_format, delimiter=table_delimiter, header=table_header, comments="")
    path_npy = ".".join(path.split(".")[:-1]) + ".npy"
    np.save(path_npy, table)


# function for creating list of format strings and a header for subsequent saving of tables into text files
def get_fmt_and_header(column_names, all_column_groups, all_data_types, delimiter="\t", precision=2, column_width=10):
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
    out = [value for value, flag in zip(iterable_to_subset, validity_mask) if flag]
    if return_array:
        return np.array(out)
    else:
        return out


# %%


## in the future, improve this so it can handle polygons etc
def check_block_for_data_and_cal(block_extents, stack_boundaries, calibration_boundaries):

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


# # %%

# def continue_with_next_block(current_block_index,block_finished_flag,block_order):
#      # swap the flag
#     block_finished_flag[current_block_index] = True
#     # remove current par block from list
#     block_order = block_order[block_order != current_block_index]
#     # select next par block
#     if (len(block_order) > 0):
#         # if there is at least one left, just take the next closest to CALdata
#         current_block_index = block_order[0]
#     else:
#         current_block_index = -1
#     return current_block_index,block_finished_flag,block_order

# %%
def compute_block_processing_order(
    block_corner_coordinates_east,
    block_corner_coordinates_north,
    block_size_east,
    block_size_north,
    calibration_area_coordinates,
    stack_data_coordinates,
):
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
    else:
        out_data = np.copy(in_data)
    return out_data


# %% function for calculating a given statistic for polygons of arbitrary shape
def stats_on_polygons(data, pixel_axis_east, pixel_axis_north, reference_polygons, method):
    # calculates statistic in method on CAL data polygons

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
    polygon_means_vec = np.zeros(len(reference_polygons))

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

        raster_support_driver = None
        poly_layer_support_driver = None

    return polygon_means_vec


# %% function for calculating a given statistic on all samples on a grid and all polygons
def stats_on_all_samples(
    data,
    pixel_axis_east,
    pixel_axis_north,
    sample_axis_east,
    sample_size_east,
    sample_axis_north,
    sample_size_north,
    polygons,
    method,
):
    stats = mean_on_rois(
        data,
        pixel_axis_east,
        pixel_axis_north,
        sample_axis_east,
        sample_size_east,
        sample_axis_north,
        sample_size_north,
        method,
    )
    if polygons:
        stats_polygons = stats_on_polygons(data, pixel_axis_east, pixel_axis_north, polygons, method)
        stats = np.concatenate((stats, stats_polygons))
    return stats


# # %%
# # function for creating list of format strings and a header for subsequent saving of tables into text files
# def get_fmt_and_header(column_names,all_column_groups,all_data_types,delimiter='\t',precision=2,column_width=10):
#     out_format = []
#     out_header = []
#     for curr_column in column_names:
#         for curr_column_group,curr_data_type in zip(all_column_groups,all_data_types):
#             if curr_column in curr_column_group:
#                 if curr_data_type == 'd':
#                     out_format.append('%s%dd' % (r'%',column_width))
#                 elif curr_data_type == 'f':
#                     out_format.append('%s%d.%df' % (r'%',column_width,precision))
#                 break
#         out_header.append('%s%ds' % (r'%',column_width) % curr_column)
#     return out_format,delimiter.join(out_header)

# #%% swap variable names in formulas to slices of an array
# def swap_names_and_merge_formula(original_formulas,observable_names,parameter_names,new_table_name,use_observable_if_repeated_and_available = True):
#     original_variable_names = observable_names + parameter_names
#     unique_variable_names,name_counts = np.unique(np.array(original_variable_names),return_counts = True)
#     new_formula = '0'
#     for current_formula in original_formulas:
#         for unique_variable_name,name_count in zip(unique_variable_names,name_counts):
#             if name_count==1:
#                 position_in_variable_names_vector = np.where(np.array(original_variable_names)==unique_variable_name)[0][0]
#             elif name_count==2:
#                 position_in_variable_names_vector = np.where(np.array(original_variable_names)==unique_variable_name)[0][np.int32(~use_observable_if_repeated_and_available)]
#             current_formula = current_formula.replace(unique_variable_name,new_table_name + ('[:,%d]' % (position_in_variable_names_vector)))
#         new_formula = new_formula + ' + (%s)**2' % (current_formula)
#     return new_formula.replace('0 + ','')
# %%
