# IMPORT NECESSARY MODULES
import logging
import numpy as np
import scipy as sp
from shapely.geometry import MultiPoint
from scipy.interpolate import interp2d
from osgeo import gdal
from gdalconst import GA_ReadOnly
from equi7grid.equi7grid import Equi7Grid
from biomassL2.utility_orchestrator import choose_equi7_sampling


def initialize_inversion_parameters(eq7_sampling, geographic_grid_sampling, geographic_boundaries, proc_confAGB):

    e7g_intermediate = Equi7Grid(eq7_sampling)

    # input pixel spacing (assumes last part of the path is the resolution
    # as is default in the simulation script)
    dx_in = eq7_sampling  # int(indir.split('_')[-1][:-1])

    # N0: upper left corner for output map (UTM 32S)
    sub_grid_string, x_upper_left, y_upper_left = e7g_intermediate.lonlat2xy(
        geographic_boundaries.lon_min, geographic_boundaries.lat_max
    )
    N0 = y_upper_left
    E0 = x_upper_left

    # lower right corner for output map (UTM 32S)
    _, x_lower_right, y_lower_right = e7g_intermediate.lonlat2xy(
        geographic_boundaries.lon_max, geographic_boundaries.lat_min
    )
    NN = y_lower_right
    EN = x_lower_right

    # ground pixel spacing for output map (m)
    if geographic_grid_sampling > (proc_confAGB.product_resolution / 2):
        geographic_grid_sampling = proc_confAGB.product_resolution / 2
        logging.warning(
            'Input "geographic_grid_sampling" cannot be greater than configuration "product_resolution"/2, setting to {}'.format(
                geographic_grid_sampling
            )
        )
    else:
        geographic_grid_sampling = geographic_grid_sampling

    geographic_grid_sampling = choose_equi7_sampling(proc_confAGB.product_resolution, geographic_grid_sampling)

    dE = np.abs(eq7_sampling)
    dN = -np.abs(eq7_sampling)

    # roi spacing for output map (m)
    dE_roi = np.abs(proc_confAGB.distance_sampling_area)
    dN_roi = -np.abs(proc_confAGB.distance_sampling_area)

    # define ROI size
    wE_roi = np.abs(proc_confAGB.product_resolution)  # dx_in # (cannot be finer than input pixel size)
    wN_roi = -np.abs(proc_confAGB.product_resolution)  # -dx_in # (cannot be finer than input pixel size)

    # parameter block spacing (m)
    dE_par = np.abs(proc_confAGB.distance_parameter_block)
    dN_par = -np.abs(proc_confAGB.distance_parameter_block)

    # parameter block width (m)
    wE_par = np.abs(
        proc_confAGB.parameter_block_size
    )  # dE_par*2 # 3/2 times spacing means that each tile shares 33.3% with its neighbour
    wN_par = -np.abs(
        proc_confAGB.parameter_block_size
    )  # dN_par*2 # 3/2 times spacing means that each tile shares 33.3% with its neighbour

    # parameter variability across acquisitions
    # columns are polarisation id, acquisition id, global cycle id, phi angle, swath id, subswath id, azimuth id
    # determine how to vary each parameter, 1 means that parameter is varied in this direction, 0 means its constant for all
    var_pol = proc_confAGB.model_parameters.changes_across_pol
    var_stack = proc_confAGB.model_parameters.changes_across_stacks
    var_gc = proc_confAGB.model_parameters.changes_across_global_cycle
    var_h = proc_confAGB.model_parameters.changes_across_heading
    var_sw = proc_confAGB.model_parameters.changes_across_swath
    var_subsw = proc_confAGB.model_parameters.changes_across_subswath
    var_az = proc_confAGB.model_parameters.changes_across_azimuth

    all_vars = [
        np.array(
            [
                var_pol.agb_model_scaling,
                var_stack.agb_model_scaling,
                var_gc.agb_model_scaling,
                var_h.agb_model_scaling,
                var_sw.agb_model_scaling,
                var_subsw.agb_model_scaling,
                var_az.agb_model_scaling,
            ]
        ),
        np.array(
            [
                var_pol.agb_model_exponent,
                var_stack.agb_model_exponent,
                var_gc.agb_model_exponent,
                var_h.agb_model_exponent,
                var_sw.agb_model_exponent,
                var_subsw.agb_model_exponent,
                var_az.agb_model_exponent,
            ]
        ),
        np.array(
            [
                var_pol.agb_cosine_exponent,
                var_stack.agb_cosine_exponent,
                var_gc.agb_cosine_exponent,
                var_h.agb_cosine_exponent,
                var_sw.agb_cosine_exponent,
                var_subsw.agb_cosine_exponent,
                var_az.agb_cosine_exponent,
            ]
        ),
    ]

    # parameter ranges
    l_lims = proc_confAGB.model_parameters.agb_model_scaling_limits
    a_lims = proc_confAGB.model_parameters.agb_model_exponent_limits
    n_lims = proc_confAGB.model_parameters.agb_cosine_exponent_limits
    W_lims = proc_confAGB.model_parameters.agb_estimation_limits
    w_lims = 10 * np.log10(W_lims)

    # number of tests to run (i.e., subsets of cals)
    N_tests = int(proc_confAGB.number_of_tests)

    # adjust ROI spacings to be intere multiples of output map ground pixel spacing
    dE_roi_old = dE_roi
    dE_roi = round(dE_roi / dE) * dE
    dN_roi = -dE_roi
    logging.info(
        'ROI ground pixel spacing (east and north) adjusted from {} to {} [m] in order to be an intere multiple of output map ground pixel spacing '.format(
            dE_roi_old, dE_roi
        )
    )

    # adjust PAR spacings to be intere multiples of ROI ground pixel spacing
    dE_par_old = dE_par
    dE_par = round(dE_par / dE_roi) * dE_roi
    dN_par = -dE_par
    logging.info(
        'PAR ground pixel spacing (east and north) adjusted from {} to {} [m] in order to be an intere multiple of ROI pixel spacing '.format(
            dE_par_old, dE_par
        )
    )

    return (
        dx_in,
        N0,
        E0,
        NN,
        EN,
        dE,
        dN,
        dE_roi,
        dN_roi,
        wE_roi,
        wN_roi,
        dE_par,
        dN_par,
        wE_par,
        wN_par,
        all_vars,
        l_lims,
        a_lims,
        n_lims,
        w_lims,
        N_tests,
        geographic_grid_sampling,
        sub_grid_string,
    )


def compute_processing_blocs_order(lut_cal, EE_par_coordinates, wE_par, NN_par_coordinates, wN_par):

    # compute cal areas and find the block in which the largest calibration area is
    num_cals = lut_cal.shape[0]
    cal_areas = np.zeros(num_cals)
    for idx_cal in np.arange(num_cals):
        east_length = lut_cal[idx_cal, 1] - lut_cal[idx_cal, 0]
        north_length = lut_cal[idx_cal, 3] - lut_cal[idx_cal, 2]
        cal_areas[idx_cal] = east_length * north_length
    largest_cal_coords = lut_cal[np.argsort(cal_areas)[-1], :-1]

    intersection_areas = compute_intersection_area(
        EE_par_coordinates,
        EE_par_coordinates + wE_par,
        NN_par_coordinates + wN_par,
        NN_par_coordinates,
        largest_cal_coords[0],
        largest_cal_coords[1],
        largest_cal_coords[2],
        largest_cal_coords[3],
    )

    par_block_index_curr = np.argsort(intersection_areas)[-1]

    # sort blocks by distance from currentblock, extract only the tiles which cover the current polygon and the entire studied area
    block_order = np.argsort(
        np.sqrt(
            (EE_par_coordinates - EE_par_coordinates[par_block_index_curr]) ** 2
            + (NN_par_coordinates - NN_par_coordinates[par_block_index_curr]) ** 2
        )
    )  # tile_order

    return par_block_index_curr, block_order


# Simple table look-up, without interpolation or extrapolation, just exact matching
def tableLookupInt(x, y, xx):
    yy = -np.int32(np.ones(xx[:, 0].shape))
    for irow in np.arange(x.shape[0]):
        yy[(xx[:, 0] == x[irow, 0]) & (xx[:, 1] == x[irow, 1])] = y[irow]
    return yy


# forward and inverse transforms for parameters
def forward(x):
    return 10 * np.log10(x)


def inverse(x):
    return 10 ** (0.1 * x)


# function for regularizing indices
def regularizeIndices(a, b, off):
    uni = np.unique(np.concatenate((a, b)))
    lut = np.column_stack((uni, off + np.arange(len(uni))))
    f = sp.interpolate.interp1d(
        np.concatenate((-np.ones(1), lut[:, 0])), np.concatenate((-np.ones(1), lut[:, 1])), kind='nearest'
    )
    return np.int32(f(a)), np.int32(f(b)), np.int32(lut)


def tableLookupFloat(x, y, xx):
    yy = np.nan * np.ones(xx[:, 0].shape)
    for irow in np.arange(x.shape[0]):
        yy[(xx[:, 0] == x[irow, 0]) & (xx[:, 1] == x[irow, 1])] = y[irow]
    return yy


def check_intersection(
    A_east_min, A_east_max, A_north_min, A_north_max, B_east_min, B_east_max, B_north_min, B_north_max
):

    if (
        (A_east_min >= A_east_max).any()
        or (A_north_min >= A_north_max).any()
        or (B_east_min >= B_east_max).any()
        or (B_north_min >= B_north_max).any()
    ):
        raise ValueError(
            'Wrong inputs, at least one of the "min" variable is greater than (or equal to) the "max" one.'
        )

    A_intersects_B = (
        (A_east_min <= B_east_max)
        & (A_east_max >= B_east_min)
        & (A_north_min <= B_north_max)
        & (A_north_max >= B_north_min)
    )

    return A_intersects_B


def compute_intersection_area(
    A_east_min_vec, A_east_max_vec, A_north_min_vec, A_north_max_vec, B_east_min, B_east_max, B_north_min, B_north_max
):

    B_coords = tuple(
        zip([B_north_min, B_north_min, B_north_max, B_north_max], [B_east_min, B_east_max, B_east_min, B_east_max])
    )
    B_polygon_obj = MultiPoint(B_coords).convex_hull

    num_areas = A_east_min_vec.shape[0]
    intersection_areas = np.zeros(num_areas)
    for area_idx in np.arange(num_areas):
        A_coords = tuple(
            zip(
                [
                    A_north_min_vec[area_idx],
                    A_north_min_vec[area_idx],
                    A_north_max_vec[area_idx],
                    A_north_max_vec[area_idx],
                ],
                [
                    A_east_min_vec[area_idx],
                    A_east_max_vec[area_idx],
                    A_east_min_vec[area_idx],
                    A_east_max_vec[area_idx],
                ],
            )
        )
        A_polygon_obj = MultiPoint(A_coords).convex_hull

        intersection_areas[area_idx] = B_polygon_obj.intersection(A_polygon_obj).area

    return intersection_areas


def interp2d_wrapper(data_path, band_index_to_read, east_out_axis, north_out_axis, fill_value):

    data_driver = gdal.Open(data_path, GA_ReadOnly)
    data_in = data_driver.GetRasterBand(band_index_to_read).ReadAsArray()
    geotransform = data_driver.GetGeoTransform()
    east_axis_in = geotransform[0] + geotransform[1] * np.arange(data_driver.RasterXSize)
    north_axis_in = geotransform[3] + geotransform[5] * np.arange(data_driver.RasterYSize)
    data_driver = None

    # interpolation function
    min_input = np.nanmin(data_in)
    nan_mask = np.isnan(data_in)
    data_in[nan_mask] = -9999

    interp2d_fun = interp2d(east_axis_in, north_axis_in[::-1], data_in[::-1, :], fill_value=fill_value)
    interp2d_fun_mask = interp2d(east_axis_in, north_axis_in[::-1], nan_mask[::-1, :], fill_value=float(0))

    # interpolate
    data_out = interp2d_fun(east_out_axis, north_out_axis[::-1])[::-1, :]
    mask_out = np.ceil(interp2d_fun_mask(east_out_axis, north_out_axis[::-1]))[::-1, :]

    mask_out = np.logical_and(mask_out, np.isnan(data_out))
    data_out[mask_out.astype(int)] = -9999

    data_out[data_out < min_input] = np.nan
    data_out[mask_out.astype(int)] = np.nan

    return data_out


def merge_agb_intermediate(a_data, b_data, method='nan_mean'):

    if method == 'nan_mean':
        a_data = np.nanmean(np.dstack((a_data, b_data)), axis=2)

    elif method == 'overlap_mean':
        # mean all the stacks
        a_data_is_present = np.logical_not(np.isnan(a_data))
        b_data_is_present = np.logical_not(np.isnan(b_data))

        # where both are present, perform the mean
        a_data[np.logical_and(b_data_is_present, a_data_is_present)] = np.mean(
            a_data[np.logical_and(b_data_is_present, a_data_is_present)],
            b_data[np.logical_and(b_data_is_present, a_data_is_present)],
        )

        # where only b is present perform the sum
        a_data[np.logical_and(b_data_is_present, np.logical_not(a_data_is_present))] = np.sum(
            a_data[np.logical_and(b_data_is_present, np.logical_not(a_data_is_present))],
            b_data[np.logical_and(b_data_is_present, np.logical_not(a_data_is_present))],
        )

    return a_data


def mean_on_rois(stack_data_interp, EE_pixel_mesh, NN_pixel_mesh, EE_roi_axis, dE_roi, NN_roi_axis, dN_roi, method):

    # Mean on ROIs
    number_of_rois = len(EE_roi_axis) * len(NN_roi_axis)
    data_roi_means_vec = np.zeros(number_of_rois)

    roi_counter = 0
    for north_roi_max in NN_roi_axis:
        for east_roi_min in EE_roi_axis:
            roi_counter = roi_counter + 1

            east_roi_max = east_roi_min + dE_roi
            north_roi_min = north_roi_max + dN_roi

            # here use < and > instead of <= and >= because the pixel meshes has been defined in the centre of the pixels
            EE_axis_valid = np.logical_and(EE_pixel_mesh > east_roi_min, EE_pixel_mesh < east_roi_max)
            NN_axis_valid = np.logical_and(NN_pixel_mesh > north_roi_min, NN_pixel_mesh < north_roi_max)
            curr_roi_indexes_pixels = np.logical_and(EE_axis_valid, NN_axis_valid)

            if method == 'mean':
                data_roi_means_vec[roi_counter - 1] = np.mean(stack_data_interp[curr_roi_indexes_pixels])
            elif method == 'nan_mean':
                data_roi_means_vec[roi_counter - 1] = np.nanmean(stack_data_interp[curr_roi_indexes_pixels])

    return data_roi_means_vec


def get_projection_from_path(file_name):

    driver = gdal.Open(file_name, GA_ReadOnly)
    projection = driver.GetProjection()

    return projection
