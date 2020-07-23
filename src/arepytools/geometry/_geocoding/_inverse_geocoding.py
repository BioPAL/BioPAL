# SPDX-FileCopyrightText: Aresys S.r.l. <info@aresys.it>
# SPDX-License-Identifier: MIT

"""
Inverse geocoding module
------------------------
"""

import copy
import numpy as np
import arepytools.constants as cst
from arepytools import _utils

_VELOCITY_SCENE = np.asarray([0, 0, 0])

_NEWTON_TOLERANCE = 1.e-8
_MAX_ITERATIONS = 8
_ZERO_TOLERANCE = 1.e-8


def inverse_geocoding_monostatic(gso, point: np.ndarray, frequency_doppler_centroid: float,
                                 wavelength: float) -> (list, list):
    """Perform monostatic inverse geocoding

    :param gso: general sar orbit
    :param point: a point (xyz) coordinate
    :param frequency_doppler_centroid: doppler centroid
    :param wavelength:  wavelength

    More than one solution is possible: that's why lists are returned.

    :return: azimuth_times, range_times. Two lists with the solutions.
    """
    _check_inverse_geocoding_input(point, frequency_doppler_centroid, wavelength)

    interval_ids = _inverse_geocoding_monostatic_init(gso, point.reshape((3, 1)), frequency_doppler_centroid,
                                                      wavelength)

    azimuth_times, range_times = list(), list()
    for interval_index in interval_ids:
        az_time, rg_time = _inverse_geocoding_monostatic_core(gso, point, frequency_doppler_centroid, wavelength,
                                                              interval_index)
        azimuth_times.append(az_time)
        range_times.append(rg_time)
    return azimuth_times, range_times


def inverse_geocoding_bistatic(gso_rx, point: np.ndarray, frequency_doppler_centroid: float,
                               wavelength: float, gso_tx) -> (list, list):
    """Perform bistatic inverse geocoding

    :param gso_rx: general sar orbit (rx)
    :param point: a point (xyz) coordinate
    :param frequency_doppler_centroid: doppler centroid
    :param wavelength:  wavelength
    :param gso_tx: general sar orbit (tx)

    More than one solution is possible: that's why lists are returned.

    :return: azimuth_times, range_times. Two lists with the solutions.
    """

    azimuth_axis_rx = gso_rx.time_axis_array
    azimuth_axis_tx = gso_tx.time_axis_array

    try:
        dt_rx = gso_rx.dt
        dt_tx = gso_tx.dt
        dt = min(dt_tx, dt_rx)
    except RuntimeError as exc:
        raise RuntimeError("Cannot perform bistatic inverse geocoding for orbit with non-regular axis") from exc

    if dt_rx <= 0 or dt_tx <= 0:
        raise RuntimeError("Cannot perform bistatic inverse geocoding for orbit with decreasing axis")

    t_min, t_max = max(azimuth_axis_rx[0], azimuth_axis_tx[0]), min(azimuth_axis_rx[-1], azimuth_axis_tx[-1])

    if t_max <= t_min:
        raise RuntimeError('Cannot perform bistatic inverse geocoding: orbits not overlapped')

    t_ref = t_min + (t_max - t_min) / 2
    t_min = t_min - t_ref
    t_max = t_max - t_ref

    t_rx_cur = 0.
    t_tx_cur = 0.
    t_rx_old = t_max
    t_tx_old = t_max

    # t_tx_cur initialization
    position_rx = gso_rx.get_position(t_ref + t_rx_cur).squeeze()
    position_tx = gso_tx.get_position(t_ref + t_tx_cur).squeeze()
    delay_estimate = (np.linalg.norm(position_tx - point) + np.linalg.norm(position_rx - point)) / cst.LIGHT_SPEED
    t_tx_cur = t_rx_cur - delay_estimate

    interval_index_rx = -1
    interval_index_tx = -1

    iteration_counter = 0
    increment_threshold = 1.e-8
    max_iterations = 8
    relative_threshold = 2.
    increment_old = 1e+20
    increment = 1e+19

    def condition(ctr: int, incr: float, incr_old: float):
        return (ctr < max_iterations) and (incr > increment_threshold) and (incr_old / incr > relative_threshold)

    while condition(iteration_counter, increment, increment_old):
        if (t_rx_cur < t_min) or (t_rx_cur > t_max) or (t_tx_cur < t_min) or (t_tx_cur > t_max):
            raise RuntimeError('Cannot perform inverse geocoding: iterations outside the time axis')

        if (abs(t_rx_old - t_rx_cur) >= dt / 2.) or (interval_index_rx < 0):
            interval_index_rx = gso_rx.interpolator.get_interval_index(t_ref + t_rx_cur)

        position_rx = gso_rx.get_position(t_ref + t_rx_cur, interval_indexes=interval_index_rx).squeeze()
        velocity_rx = gso_rx.get_velocity(t_ref + t_rx_cur, interval_indexes=interval_index_rx).squeeze()
        acceleration_rx = gso_rx.get_acceleration(t_ref + t_rx_cur, interval_indexes=interval_index_rx).squeeze()

        if (abs(t_tx_old - t_tx_cur) >= dt / 2.) or (interval_index_tx < 0):
            interval_index_tx = gso_tx.interpolator.get_interval_index(t_ref + t_tx_cur)

        position_tx = gso_tx.get_position(t_ref + t_tx_cur, interval_indexes=interval_index_tx).squeeze()
        velocity_tx = gso_tx.get_velocity(t_ref + t_tx_cur, interval_indexes=interval_index_tx).squeeze()
        acceleration_tx = gso_tx.get_acceleration(t_ref + t_tx_cur, interval_indexes=interval_index_tx).squeeze()

        norm_vel_rx_square = np.dot(velocity_rx, velocity_rx.T)
        norm_vel_tx_square = np.dot(velocity_tx, velocity_tx.T)

        slant_range_rx = np.linalg.norm(position_rx - point)
        slant_range_tx = np.linalg.norm(position_tx - point)
        slant_range = (t_rx_cur - t_tx_cur) * cst.LIGHT_SPEED

        range_dot_vel_rx = np.dot(point - position_rx, velocity_rx)
        range_dot_vel_tx = np.dot(point - position_tx, velocity_tx)

        distance_equation_residual = slant_range_rx + slant_range_tx - slant_range
        doppler_equation_residual = -range_dot_vel_tx * slant_range_rx - range_dot_vel_rx * slant_range_tx
        doppler_equation_residual += wavelength * frequency_doppler_centroid * slant_range_rx * slant_range_tx

        df1_dt_tx = range_dot_vel_tx / slant_range_tx + cst.LIGHT_SPEED
        df1_dt_rx = range_dot_vel_rx / slant_range_rx - cst.LIGHT_SPEED
        df2_dt_tx = slant_range_rx * (norm_vel_tx_square - np.dot(point - position_tx, acceleration_tx))
        df2_dt_tx += range_dot_vel_tx / slant_range_tx * (
            range_dot_vel_rx + wavelength * frequency_doppler_centroid * slant_range_rx)
        df2_dt_rx = slant_range_tx * (norm_vel_rx_square - np.dot(point - position_rx, acceleration_rx))
        df2_dt_rx += range_dot_vel_rx / slant_range_rx * (
            range_dot_vel_tx + wavelength * frequency_doppler_centroid * slant_range_tx)

        jac_det = df1_dt_tx * df2_dt_rx - df1_dt_rx * df2_dt_tx
        inv_jac = 1. / jac_det * np.asarray([[df2_dt_rx, -df1_dt_rx], [-df2_dt_tx, df1_dt_tx]])
        delta = inv_jac @ np.asarray([distance_equation_residual, doppler_equation_residual])

        t_tx_old = t_tx_cur
        t_rx_old = t_rx_cur

        t_tx_cur = t_tx_cur - delta[0]
        t_rx_cur = t_rx_cur - delta[1]

        increment_old = increment
        increment = np.dot(delta, delta)

        iteration_counter += 1

    # lists for compatibility with monostatic version.
    t_azimuth = [t_ref + t_rx_cur]
    t_range = [t_rx_cur - t_tx_cur]
    return t_azimuth, t_range


def _inverse_geocoding_monostatic_init(gso, point: np.ndarray, frequency_doppler_centroid: float, wavelength: float):
    """Function to compute initialization of inverse geocoding

    The azimuth time we are looking for is such that the doppler equation is equal to zero.
    In this function we look for the time interval where such equation crosses zero. We compute the corresponding
    polynomial index that will be used subsequently for the precise computation of the aximuth time.
    This function evaluates the doppler equation on the azimuth time axis of the gso. It checks where it changes sign
    (from plus to minus). For each of those time positions the index of polynomial to use for interpolating in that
    neighbourhood is computed and returned. If no crossing points are found, either the first or the last polynomials
    ids are returned depending on the time instants that has smaller values of the doppler equation.

    :param gso: general sar orbit
    :param point: xyz point (ecef)
    :param frequency_doppler_centroid: a float with the doppler centroid
    :param wavelength: the wavelenth

    :return: (a list of) index of the intervals to be used for the inverse geocoding.
    """
    doppler_centroid_equation = gso.evaluate_doppler_equation(point, frequency_doppler_centroid, wavelength)

    zero_crossing_indexes = list()
    for k in range(1, doppler_centroid_equation.size):
        if doppler_centroid_equation[k] * doppler_centroid_equation[k - 1] <= 0 and doppler_centroid_equation[k] < \
                doppler_centroid_equation[k - 1]:
            zero_crossing_indexes.append(k)

    interval_index = list()
    if zero_crossing_indexes:
        azimuth_time = gso.get_interpolated_time_axis(np.asarray(zero_crossing_indexes) - 0.5)
        interval_index = gso.interpolator.get_interval_index(azimuth_time).tolist()
    else:
        if abs(doppler_centroid_equation[0]) < abs(doppler_centroid_equation[-1]):
            interval_index.append(0)
        else:
            interval_index.append(gso.interpolator.num_polynomials - 1)
    return interval_index


def _inverse_geocoding_monostatic_core(gso, point: np.ndarray, frequency_doppler_centroid: float, wavelength: float,
                                       interval_index: int):
    """Core algorithm for precise inverse geocoding

    This function finds one zero of the doppler equation as a function of the azimuth time.
    The search of the zero is restricted to the interval associated to the polynomial indexed by polynomial_index.
    Newton iterations are performed to find the exact azimuth time, once that the sensor position is known the
    computation of the range time is trivial.

    :param gso: general sar orbit
    :param point: xyz point (ecef)
    :param frequency_doppler_centroid: a float with the doppler centroid
    :param wavelength: the wavelenth
    :param interval_index: interval where to search for the solution.

    :return: (azimuth_time, range_time)

    See also _inverse_geocoding_monostatic_init
    """
    counter = 0

    polynomial_index = gso.interpolator.get_polynomial_index(interval_indexes=interval_index)[0]
    t_min = gso.time_axis_array[polynomial_index]
    t_max = gso.time_axis_array[polynomial_index + gso.interpolator.get_polynomial_order()]

    time_old = t_min
    current_azimuth_time = t_min + (t_max - t_min) / 2
    slant_range = 0
    while counter < _MAX_ITERATIONS and abs(current_azimuth_time - time_old) > _NEWTON_TOLERANCE:
        if current_azimuth_time < t_min or current_azimuth_time > t_max:
            raise RuntimeError('Cannot perform inverse geocoding')
        time_old = copy.deepcopy(current_azimuth_time)
        sat_position = gso.get_position(current_azimuth_time, interval_indexes=interval_index)[:, 0]
        sat_velocity = gso.get_velocity(current_azimuth_time, interval_indexes=interval_index)[:, 0]
        sat_acceleration = gso.get_acceleration(current_azimuth_time, interval_indexes=interval_index)[:, 0]

        point2sat = point - sat_position
        slant_range = np.linalg.norm(point2sat)
        f = np.dot(point2sat,
                   _VELOCITY_SCENE - sat_velocity) + wavelength * frequency_doppler_centroid / 2. * slant_range
        df = - np.dot(sat_velocity, _VELOCITY_SCENE - sat_velocity) - np.dot(sat_acceleration, point2sat) \
             - wavelength * frequency_doppler_centroid / 2 * np.dot(sat_velocity, point2sat) / slant_range
        current_azimuth_time -= f / df

        counter = counter + 1

    if counter == 8:
        raise RuntimeError("Exceeded max number of iterations during inverse geocoding")

    return current_azimuth_time, slant_range * 2 / cst.LIGHT_SPEED


def _check_inverse_geocoding_input(point, frequency_doppler_centroid, wavelength):
    """Checks that the input parameters for the inverse geocoding have proper types and sizes

    :param point: should be any 3-sized np.ndarray
    :param frequency_doppler_centroid: should be a float
    :param wavelength: should be a float
    """

    # checking types
    _utils.check_type(point, np.ndarray, "point")
    _utils.check_type(frequency_doppler_centroid, float, "frequency_doppler_centroid")
    _utils.check_type(wavelength, float, "wavelength")

    # size check
    if point.size != 3:
        RuntimeError("Point should be a np.ndarray of size 3")
