# SPDX-FileCopyrightText: Aresys S.r.l. <info@aresys.it>
# SPDX-License-Identifier: MIT

"""
Direct geocoding module
-----------------------
"""

import numpy as np

from arepytools.geometry import conversions as conv
from arepytools.geometry.wgs84 import WGS84
import arepytools.constants as cst
from arepytools.geometry._geocoding import _newton
from arepytools import _utils

_LOOK_SIGN = {"RIGHT": 1, "LEFT": -1}


def _mid_range_distance(range_axis: np.ndarray):
    return range_axis[(range_axis.size + 1) // 2 - 1] * cst.LIGHT_SPEED / 2


def direct_geocoding_monostatic(sensor_position: np.ndarray, sensor_velocity: np.ndarray,
                                range_times: np.ndarray, look_direction: str, geodetic_altitude: float,
                                frequency_doppler_centroid: np.ndarray, wavelength: float,
                                initial_guess: np.ndarray = None) -> np.ndarray:
    """Perform monostatic direct geocoding

    :param sensor_position: position of the sensor
    :param sensor_velocity: velocity of the sensor
    :param range_times: range axis (N, 1) vector
    :param look_direction: either RIGHT or LEFT
    :param geodetic_altitude: the altitude over wgs84
    :param frequency_doppler_centroid: array with frequency doppler centroid values (N,1) vector
    :param wavelength: the wavelength
    :param initial_guess: initial guess for newton iterations. If not provided a guess will be computed [optional]

    :return: a matrix 3xN with the xyz coordinate of the points
    """

    # Input validation
    _utils.check_type(wavelength, float)
    _utils.check_type(geodetic_altitude, float)
    sensor_position = _utils.input_data_to_numpy_array_with_checks(sensor_position, dtype=float, shape=(3,))
    sensor_velocity = _utils.input_data_to_numpy_array_with_checks(sensor_velocity, dtype=float, shape=(3,))

    # Optional computation of the initial guess
    if initial_guess is None:
        initial_guess = _direct_geocoding_monostatic_init(sensor_position, sensor_velocity,
                                                          _mid_range_distance(range_times),
                                                          _LOOK_SIGN[look_direction])

    # Geocoding
    output = np.zeros((3, range_times.size))
    for i_rg, (rg_time, f_dc) in enumerate(zip(range_times, frequency_doppler_centroid)):
        output[:, i_rg] = _direct_geocoding_monostatic_core(initial_guess, sensor_position, sensor_velocity, rg_time,
                                                            wavelength, f_dc, geodetic_altitude)
        initial_guess = output[:, i_rg].copy()
    return output


def direct_geocoding_bistatic(sensor_position_rx: np.ndarray, sensor_velocity_rx: np.ndarray,
                              sensor_positions_tx: np.ndarray, sensor_velocities_tx: np.ndarray,
                              range_times: np.ndarray,
                              look_direction: str, geodetic_altitude: float,
                              frequency_doppler_centroid: np.ndarray, wavelength: float,
                              initial_guess: np.ndarray = None) -> np.ndarray:
    """Perform bistatic direct geocoding

    :param sensor_position_rx: position of the sensor rx
    :param sensor_velocity_rx: velocity of the sensor rx
    :param sensor_positions_tx: positions of the sensor tx (one for each range time (3, N))
    :param sensor_velocities_tx: velocities of the sensor tx (one for each range time (3, N))
    :param range_times: range axis (N, 1) vector
    :param look_direction: either RIGHT or LEFT
    :param geodetic_altitude: the altitude over wgs84
    :param frequency_doppler_centroid: array with frequency doppler centroid values (N,1) vector
    :param wavelength: the wavelength
    :param initial_guess: initial guess for newton iterations. If not provided a guess will be computed [optional]

    :return: a matrix 3xN with the xyz coordinate of the points
    """
    # input validation
    _utils.check_type(wavelength, float)
    _utils.check_type(geodetic_altitude, float)
    sensor_position_rx = _utils.input_data_to_numpy_array_with_checks(sensor_position_rx, dtype=float, shape=(3,))
    sensor_velocity_rx = _utils.input_data_to_numpy_array_with_checks(sensor_velocity_rx, dtype=float, shape=(3,))
    sensor_positions_tx = _utils.input_data_to_numpy_array_with_checks(sensor_positions_tx, dtype=float,
                                                                       first_axis_size=3)
    sensor_velocities_tx = _utils.input_data_to_numpy_array_with_checks(sensor_velocities_tx, dtype=float,
                                                                        first_axis_size=3)

    # Optional initial guess
    if initial_guess is None:
        initial_guess = _direct_geocoding_monostatic_init(sensor_position_rx, sensor_velocity_rx,
                                                          _mid_range_distance(range_times),
                                                          _LOOK_SIGN[look_direction])

    output = np.zeros((3, range_times.size))
    for i_rg, (rg_time, pos_tx, vel_tx, f_dc) in enumerate(
            zip(range_times, sensor_positions_tx.T, sensor_velocities_tx.T, frequency_doppler_centroid)):
        output[:, i_rg] = _direct_geocoding_bistatic_core(initial_guess, sensor_position_rx, sensor_velocity_rx, pos_tx,
                                                          vel_tx, rg_time, wavelength, f_dc,
                                                          geodetic_altitude)
        initial_guess = output[:, i_rg].copy()
    return output


def _direct_geocoding_monostatic_init(sat_position, sat_velocity, range_distance, look_sign):
    satellite_distance_from_center = np.linalg.norm(sat_position)  # Satellite Position Norm
    llh_sat = conv.xyz2llh(sat_position)
    earth_radius = np.linalg.norm(conv.llh2xyz(np.asarray([llh_sat[0], llh_sat[1], 0], dtype=float)))

    # Check earth radius vs range compatibility
    if range_distance < satellite_distance_from_center - earth_radius:
        raise RuntimeError("Cannot find initial guess for direct geocoding")

    ux = sat_position / satellite_distance_from_center
    uy = np.cross(sat_position, sat_velocity)
    uy = uy / np.linalg.norm(uy)
    uz = np.cross(ux, uy)

    # x-coordinate
    x = (satellite_distance_from_center ** 2 + earth_radius ** 2 - range_distance ** 2) / (
        2 * satellite_distance_from_center)

    # Circle radius
    ro = np.sqrt(earth_radius ** 2 - x ** 2)

    # Project velocity on ref frame
    vx = np.dot(sat_velocity, ux)
    vz = np.dot(sat_velocity, uz)

    # Find first solution
    z = (satellite_distance_from_center - x) * vx / vz
    y = np.sqrt(ro ** 2 - z ** 2)

    if look_sign * y > 0:
        y = -y
    return x * ux + y * uy + z * uz


def _direct_geocoding_monostatic_core(initial_guess, sat_position, sat_velocity, range_time, wavelength,
                                      frequency_doppler_centroid, geodetic_altitude):
    d_range_square = (cst.LIGHT_SPEED * range_time / 2.) ** 2

    geoid_r_min = WGS84.semi_minor_axis + geodetic_altitude
    geoid_r_max = WGS84.semi_major_axis + geodetic_altitude

    r_ep2 = geoid_r_min ** 2
    r_ee2 = geoid_r_max ** 2

    def direct_geocoding_function(x):
        sat2point, distance_square, distance, pv_scalar = _basic_data(x, sat_position, sat_velocity)

        range_equation = distance_square - d_range_square
        grad_range_equation = - 2 * sat2point

        doppler_equation, grad_doppler_equation = _doppler_equation(wavelength, pv_scalar, distance,
                                                                    frequency_doppler_centroid, sat_velocity,
                                                                    sat2point)

        fun = [range_equation, _ellipse_equation(x, r_ee2, r_ep2), doppler_equation]
        jacobian = [[grad_range_equation[k], _der_ellipse_equation_xi(x, k, r_ee2, r_ep2), grad_doppler_equation[k]]
                    for
                    k in range(3)]

        return fun, jacobian

    return _newton.newton_for_geocoding(direct_geocoding_function, initial_guess)


def _direct_geocoding_bistatic_core(initial_guess, sat_position_rx, sat_velocity_rx, sat_position_tx,
                                    sat_velocity_tx,
                                    range_time, wavelength, frequency_doppler_centroid, geodetic_altitude):
    d_range_square = (cst.LIGHT_SPEED * range_time) ** 2  # two-way distance

    geoid_r_min = WGS84.semi_minor_axis + geodetic_altitude
    geoid_r_max = WGS84.semi_major_axis + geodetic_altitude

    r_ep2 = geoid_r_min ** 2
    r_ee2 = geoid_r_max ** 2

    def direct_geocoding_function(x):
        sat2point_rx, _, distance_rx, pv_scalar_rx = _basic_data(x, sat_position_rx,
                                                                 sat_velocity_rx)
        sat2point_tx, _, distance_tx, pv_scalar_tx = _basic_data(x, sat_position_tx,
                                                                 sat_velocity_tx)

        # range equation
        distance = distance_rx + distance_tx
        range_equation = distance ** 2 - d_range_square
        grad_range_equation = - 2 * distance * (sat2point_rx / distance_rx + sat2point_tx / distance_tx)

        # doppler equation
        doppler_equation_rx, grad_doppler_equation_rx = _doppler_equation(wavelength, pv_scalar_rx, distance_rx,
                                                                          frequency_doppler_centroid,
                                                                          sat_velocity_rx,
                                                                          sat2point_rx)
        doppler_equation_tx, grad_doppler_equation_tx = _doppler_equation(wavelength, pv_scalar_tx, distance_tx,
                                                                          frequency_doppler_centroid,
                                                                          sat_velocity_tx,
                                                                          sat2point_tx)

        doppler_equation = (doppler_equation_rx + doppler_equation_tx) / 2
        grad_doppler_equation = (grad_doppler_equation_rx + grad_doppler_equation_tx) / 2

        fun = [range_equation, _ellipse_equation(x, r_ee2, r_ep2), doppler_equation]
        jacobian = [[grad_range_equation[k], _der_ellipse_equation_xi(x, k, r_ee2, r_ep2), grad_doppler_equation[k]]
                    for
                    k in range(3)]

        return fun, jacobian

    return _newton.newton_for_geocoding(direct_geocoding_function, initial_guess)


def _ellipse_equation(x, r_ee2, r_ep2):
    return (x[0] * x[0] + x[1] * x[1]) / r_ee2 + x[2] * x[2] / r_ep2 - 1.0


def _der_ellipse_equation_xi(x, i_coord, r_ee2, r_ep2):
    r2 = r_ee2 if i_coord < 2 else r_ep2
    return 2 * x[i_coord] / r2


def _doppler_equation(wavelength, pv_scalar, distance, frequency_doppler_centroid, sat_velocity, sat2point):
    c = 2. / wavelength / distance
    doppler_equation = c * pv_scalar + frequency_doppler_centroid
    grad_doppler_equation = c * (- sat_velocity + pv_scalar * sat2point / distance ** 2)
    return doppler_equation, grad_doppler_equation


def _basic_data(x, sat_position, sat_velocity):
    sat2point = sat_position - x
    distance_square = np.dot(sat2point, sat2point)
    distance = np.sqrt(distance_square)
    pv_scalar = np.dot(sat_velocity, sat2point)
    return sat2point, distance_square, distance, pv_scalar
