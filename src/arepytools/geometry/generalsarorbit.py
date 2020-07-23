# SPDX-FileCopyrightText: Aresys S.r.l. <info@aresys.it>
# SPDX-License-Identifier: MIT

"""
General SAR orbit module
------------------------
"""

import os
from typing import Union
import numpy as np

from arepytools.math import axis as are_ax
from arepytools import _utils
from arepytools.timing.precisedatetime import PreciseDateTime
from arepytools.io.metadata import StateVectors
from arepytools.geometry._interpolator import GeometryInterpolator
from arepytools.geometry._geocoding import _direct_geocoding, _inverse_geocoding


class GeneralSarOrbit:
    """
    GeneralSarOrbit class
    """

    _MINIMUM_NUMBER_OF_DATA_POINTS = GeometryInterpolator.get_min_num_of_data_points()

    @classmethod
    def get_minimum_number_of_data_points(cls):
        """Return the required minimum number of orbit data points

        :return: required minimum number of orbit data points
        """
        return cls._MINIMUM_NUMBER_OF_DATA_POINTS

    @property
    def t0(self):
        """Origin of the time axis"""
        return self._time_axis.start

    @property
    def dt(self):
        """Time axis step (if applicable)"""
        if isinstance(self._time_axis, are_ax.RegularAxis):
            return self._time_axis.step
        raise RuntimeError("Time step is not available for orbits constructed with non-regular time axis")

    @property
    def n(self):
        """Number of orbit data points"""
        return self._time_axis.size

    @property
    def position_sv(self):
        """State vectors as 1D array of size 3N"""
        return self._state_vectors.reshape((self.n * 3,), order='F')

    @property
    def time_axis_array(self):
        """Time axis as array of time points"""
        return self._time_axis.get_array()

    @property
    def interpolator(self):
        """Geometry interpolator object"""
        return self._interpolator

    def __init__(self, time_axis: Union[np.ndarray, are_ax.Axis], state_vectors: np.ndarray):
        """Initialize the general sar orbit object using the provided input parameters

        :param time_axis: time axis of length N, as Axis or numpy array of PreciseDateTime objects.
        :param state_vectors: state vectors as 1D numpy array of size 3N
        """

        if isinstance(time_axis, np.ndarray):
            if time_axis.dtype == PreciseDateTime:
                time_axis_start = time_axis[0]
                relative_time_axis = (time_axis - time_axis_start).astype(float)
                self._time_axis = are_ax.Axis(relative_time_axis, time_axis_start)
            else:
                raise ValueError("Axis should be a vector of PreciseDateTime objects")
        else:
            self._time_axis = time_axis

        _check_init_input(self._time_axis, state_vectors)

        # state_vector are stored as (3, N) numpy array
        self._state_vectors = np.vstack((state_vectors[::3], state_vectors[1::3], state_vectors[2::3]))
        self._interpolator = GeometryInterpolator(self._time_axis, self._state_vectors)

    def get_position(self, time_points, interval_indexes=None):
        """Return the sensor positions at the specified time points

        :param time_points: 1D numpy array of length N of absolute time points
        :param interval_indexes: intervals of the time axis where the given time points are expected [optional]

        :return: 3xN numpy array of sensor positions
        """
        return self.interpolator.eval(time_points, interval_indexes)

    def get_velocity(self, time_points, interval_indexes=None):
        """Return the sensor velocities at the specified time points, evaluated using the first derivative of the
        position data

        :param time_points: 1D numpy array of length N of absolute time points
        :param interval_indexes: intervals of the time axis where the given time points are expected [optional]

        :return: 3xN numpy array of sensor velocity vectors
        """
        return self.interpolator.eval_first_derivative(time_points, interval_indexes)

    def get_acceleration(self, time_points, interval_indexes=None):
        """Return the sensor accelerations at the specified time points, evaluated using the second derivative of the
        position data

        :param time_points: 1D numpy array of length N of absolute time points
        :param interval_indexes: intervals of the time axis where the given time points are expected [optional]

        :return: 3xN numpy array of sensor acceleration vectors
        """
        return self.interpolator.eval_second_derivative(time_points, interval_indexes)

    def sat2earth(self, time_point, range_times, look_direction, geodetic_altitude=0.,
                  doppler_centroid=None, carrier_wavelength=None, orbit_tx=None, bistatic_delay=True):
        """Compute monostatic or bistatic sensor-to-earth projection

        :param time_point: absolute time point
        :param range_times: a vector (N, 1) of range times or a single range time, as numpy array or scalar
        :param look_direction: either 'LEFT' or 'RIGHT'
        :param geodetic_altitude: geodetic altitude over WGS84 [optional]
        :param doppler_centroid: doppler centroid frequency as scalar or numpy array of size (N X 1) [optional]
        :param carrier_wavelength: sensor carrier wavelength [optional]
        :param orbit_tx: orbit of the TX sensor (for bistatic case) [optional]
        :param bistatic_delay: set it to false to evaluate TX sensor position at the RX time (for bistatic case)
                               [optional]

        :return: 3xN numpy array of earth positions in xyz coordinates
        """
        range_times_checked, doppler_centroid_checked = _check_sat2earth_input(time_point, range_times,
                                                                               doppler_centroid, carrier_wavelength)

        doppler_centroid_array = np.zeros(
            range_times_checked.shape) if doppler_centroid_checked is None else doppler_centroid_checked
        carrier_wavelength = 1. if carrier_wavelength is None else carrier_wavelength

        position_rx = self.get_position(time_point).squeeze()
        velocity_rx = self.get_velocity(time_point).squeeze()
        if orbit_tx is None:
            return _direct_geocoding.direct_geocoding_monostatic(position_rx, velocity_rx, range_times_checked,
                                                                 look_direction,
                                                                 geodetic_altitude, doppler_centroid_array,
                                                                 carrier_wavelength)
        else:
            if bistatic_delay:
                tx_time_points = [time_point - range_time for range_time in range_times_checked]
            else:
                tx_time_points = [time_point for _ in range_times_checked]

            positions_tx = orbit_tx.get_position(tx_time_points)
            velocities_tx = orbit_tx.get_velocity(tx_time_points)

            return _direct_geocoding.direct_geocoding_bistatic(position_rx, velocity_rx, positions_tx, velocities_tx,
                                                               range_times_checked, look_direction, geodetic_altitude,
                                                               doppler_centroid_array, carrier_wavelength)

    def earth2sat(self, earth_point, doppler_centroid=None, carrier_wavelength=None, orbit_tx=None):
        """Compute monostatic or bistatic earth-to-sat projection

        :param earth_point: xyz coordinates of the point on earth as 3D numpy array
        :param doppler_centroid: doppler centroid frequency [optional]
        :param carrier_wavelength: sensor carrier wavelength [optional]
        :param orbit_tx: orbit of the TX sensor (for bistatic case) [optional]

        :return: sat projections as two lists of azimuth_times and range_times
        """
        if doppler_centroid is None:
            doppler_centroid = 0.0
        if carrier_wavelength is None:
            carrier_wavelength = 1.0

        _check_earth2sat_input(earth_point)

        if orbit_tx is None:
            return _inverse_geocoding.inverse_geocoding_monostatic(self, earth_point, doppler_centroid,
                                                                   carrier_wavelength)
        else:
            return _inverse_geocoding.inverse_geocoding_bistatic(self, earth_point, doppler_centroid,
                                                                 carrier_wavelength,
                                                                 orbit_tx)

    def evaluate_doppler_equation(self, earth_point: np.ndarray, doppler_centroid, carrier_wavelength):
        """Evaluate the doppler equation over the orbit time axis

        .. math::

            \\frac{2}{\\lambda} \\frac{(P_{sat}(t) - P_0) \\cdot V_{sat}(t)}{||P_{sat}(t) - P_0||} - f_{DC}

        :param earth_point: xyz coordinates of the point on earth as 3D numpy array
        :param doppler_centroid: doppler centroid frequency
        :param carrier_wavelength: sensor carrier wavelength

        :return: resulting doppler equation values in a numpy array of size N
        """
        return doppler_equation(earth_point, self._state_vectors, self.get_velocity(self._time_axis.get_array()),
                                doppler_centroid, carrier_wavelength)

    def __repr__(self):
        axis_str = str(self._time_axis)
        state_vec_str = str(self._state_vectors)

        axis_portion = "Orbit defined on azimuth axis: " + os.linesep + axis_str + os.linesep
        state_vectors_portion = "State vectors: " + os.linesep + state_vec_str + os.linesep
        return axis_portion + state_vectors_portion

    def get_interpolated_time_axis(self, interpolation_positions):
        """Return the interpolated time axis on the given interpolation positions"""
        return self._time_axis.interpolate(interpolation_positions)


def create_general_sar_orbit(state_vectors: StateVectors):
    """Create general sar orbit object from state vectors metadata

    :param state_vectors: state vectors as a StateVectors metadata object

    :return: the new GeneralSarOrbit object
    """
    time_axis = are_ax.RegularAxis((0, state_vectors.time_step, state_vectors.number_of_state_vectors),
                                   state_vectors.reference_time)
    return GeneralSarOrbit(time_axis, state_vectors.position_vector.reshape((state_vectors.position_vector.size,)))


def doppler_equation(point, sensor_position, sensor_velocity, frequency_doppler_centroid, wavelength):
    """Evaluate doppler equation"""
    point2sensor = point - sensor_position
    distance = np.linalg.norm(point2sensor, axis=0)

    def col_wise_scalar_product(matrix_a, matrix_b):
        return np.einsum('ij,ij->j', matrix_a, matrix_b)  # Einstein notation -- col wise dot product.

    return np.divide(2 / wavelength * col_wise_scalar_product(point2sensor, sensor_velocity),
                     distance) - frequency_doppler_centroid


def _check_init_input(axis, state_vectors: np.ndarray):
    if not isinstance(axis, are_ax.Axis):
        raise RuntimeError("Axis should be of type {} != {}".format(are_ax.Axis, type(axis)))
    if not isinstance(axis.start, PreciseDateTime):
        raise RuntimeError("Wrong axis start: {} != {}".format(type(axis.start), PreciseDateTime))

    _utils.check_type(state_vectors, np.ndarray, "state_vectors")

    if axis.size != state_vectors.size / 3:
        raise RuntimeError("Size of state vectors not compatible with size of time axis: {} != {}".format(
            axis.size, state_vectors.size / 3))

    if axis.size < GeneralSarOrbit.get_minimum_number_of_data_points():
        raise RuntimeError(
            "Not enough state vectors provided. {} < {}".format(axis.size,
                                                                GeneralSarOrbit.get_minimum_number_of_data_points()))

    _utils.check_dtype_of_numpy_array(state_vectors, float, name='state_vectors')
    _utils.check_shape_of_numpy_array(state_vectors, (state_vectors.size,), name='state_vectors')


def _check_sat2earth_input(azimuth_time, range_times, frequency_doppler_centroid, wavelength):
    if not isinstance(azimuth_time, PreciseDateTime):
        raise RuntimeError("Azimuth should be a single absolute time")

    range_times = _utils.input_data_to_numpy_array_with_checks(range_times, dtype=float)

    if (frequency_doppler_centroid is not None and wavelength is None) or (
            wavelength is not None and frequency_doppler_centroid is None):
        raise RuntimeError("Frequency doppler centroid and wavelength should be both specified")

    if frequency_doppler_centroid is not None:
        frequency_doppler_centroid = _utils.input_data_to_numpy_array_with_checks(frequency_doppler_centroid,
                                                                                  dtype=float, ndim=1)

    if frequency_doppler_centroid is not None and (
            frequency_doppler_centroid.ndim != 1 or frequency_doppler_centroid.shape != range_times.shape):
        if frequency_doppler_centroid.size == 1:
            frequency_doppler_centroid = np.full(range_times.shape, frequency_doppler_centroid[0])
        else:
            raise RuntimeError(
                "Frequency doppler centroid vector should have the same shape of the range times vector: "
                + "{} != {}".format(frequency_doppler_centroid.shape, range_times.shape))

    return range_times, frequency_doppler_centroid


def _check_earth2sat_input(point):
    _utils.check_type(point, np.ndarray)
    _utils.check_dtype_of_numpy_array(point, float)
    _utils.check_shape_of_numpy_array(point, shape=(3,))


def _check_numpy_array_of_mjd(vec):
    _utils.check_type(vec, np.ndarray)
    _utils.check_ndim_of_numpy_array(vec, 1)
    _utils.check_dtype_of_numpy_array(vec, PreciseDateTime)
