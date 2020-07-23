"""
General SAR attitude module
---------------------------
"""

import os
from typing import Union
from enum import Enum
import numpy as np

from arepytools import _utils
from arepytools.math import axis as are_ax
from arepytools.geometry.generalsarorbit import GeneralSarOrbit, create_general_sar_orbit
from arepytools.io.metadata import AttitudeInfo, StateVectors
from arepytools.geometry._interpolator import GeometryInterpolator


class ReferenceFrame(Enum):
    """Attitude reference frame"""
    geocentric = 'GEOCENTRIC'
    geodetic = 'GEODETIC'
    zero_doppler = 'ZERODOPPLER'


class Angles(Enum):
    """Attitude angles (Yaw, Pitch, Roll)"""
    Y = 0
    P = 1
    R = 2


class RotationOrder(Enum):
    """Yaw / Pitch / Roll rotation order"""
    ypr = 'YPR'
    yrp = 'YRP'
    pry = 'PRY'
    pyr = 'PYR'
    ryp = 'RYP'
    rpy = 'RPY'


class EnumerateRot:
    """Iterate on angles of a given a rotation order.

    Example:

    .. code-block:: python

        for angleId in EnumerateRot(RotationOrder('PRY')):
            print("{}, corresponding id: {}".format(Angles(angleId), angleId))
    """

    def __init__(self, rotation_order: RotationOrder):
        """Initialize the object to iterate on the given rotation order"""
        self._rotation_order = rotation_order.value
        self._current = 0

    def __iter__(self):
        return self

    def __next__(self):
        if self._current < 3:
            self._current += 1
            return Angles[self._rotation_order[self._current - 1]].value
        else:
            raise StopIteration()


class GeneralSarAttitude:
    """
    GeneralSarAttitude class
    """

    _MINIMUM_NUMBER_OF_DATA_POINTS = GeometryInterpolator.get_min_num_of_data_points()

    @classmethod
    def get_minimum_number_of_data_points(cls):
        """Return the required minimum number of attitude data points

        :return: required minimum number of attitude data points
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
        raise RuntimeError("Time step is not available for attitudes constructed with non-regular time axis")

    @property
    def n(self):
        """Number of attitude data points"""
        return self._time_axis.size

    @property
    def ypr_angles(self):
        """YPR angles as a 3xN numpy array"""
        return self._ypr_angles

    @property
    def rotation_order(self):
        """Yaw/Pitch/Roll rotation order"""
        return self._rotation_order

    @property
    def reference_frame(self):
        """Reference frame"""
        return self._reference_frame

    @property
    def time_axis_array(self):
        """Time axis as array of time points"""
        return self._time_axis.get_array()

    @property
    def interpolator(self):
        """Geometry interpolator object"""
        return self._interpolator

    def __init__(self, orbit: GeneralSarOrbit, time_axis: Union[np.ndarray, are_ax.Axis],
                 ypr_angles: np.ndarray, rotation_order: str, reference_frame: str):
        """Initialize the general sar attitude using the provided input parameters

        :param orbit: sensor orbit as a general sar orbit object
        :param time_axis: time axis of length N, as Axis or numpy array of PreciseDateTime objects.
        :param ypr_angles: yaw, pitch and roll angles as 3xN numpy array
        :param rotation_order: rotation order string
        :param reference_frame: reference frame string
        """

        if isinstance(time_axis, np.ndarray):
            time_axis_start = time_axis[0]
            relative_time_axis = (time_axis - time_axis_start).astype(float)
            self._time_axis = are_ax.Axis(relative_time_axis, time_axis_start)
        else:
            self._time_axis = time_axis

        _check_init_input(orbit, self._time_axis, ypr_angles, rotation_order, reference_frame)

        self._orbit = orbit
        self._rotation_order = RotationOrder(rotation_order.upper())
        self._reference_frame = ReferenceFrame(reference_frame.upper())
        self._ypr_angles = ypr_angles

        self._interpolator = GeometryInterpolator(self._time_axis, self._ypr_angles)

    def _interpolate_angles(self, time_points, angles_to_interpolate, interval_indexes):
        """Interpolate the required angle components on the given time points

        :param time_points: 1D numpy array of length N of absolute time points
        :param angles_to_interpolate: C-length vector of indexes of the angle components to interpolate
        :param interval_indexes: intervals of the time axis where the given time points are expected [optional]

        :return: CxN numpy array of interpolated angle values"""
        return self.interpolator.eval(time_points, interval_indexes, angles_to_interpolate)

    def get_yaw(self, time_points, interval_indexes=None):
        """Return the yaw angles at the specified time points

        :param time_points: 1D numpy array of length N of absolute time points
        :param interval_indexes: intervals of the time axis where the given time points are expected [optional]

        :return: 1D numpy array of N yaw angles
        """
        return self._interpolate_angles(time_points, [Angles.Y.value], interval_indexes)

    def get_pitch(self, time_points, interval_indexes=None):
        """Return the pitch angles at the specified time points

        :param time_points: 1D numpy array of length N of absolute time points
        :param interval_indexes: intervals of the time axis where the given time points are expected [optional]

        :return: 1D numpy array of N pitch angles
        """
        return self._interpolate_angles(time_points, [Angles.P.value], interval_indexes)

    def get_roll(self, time_points, interval_indexes=None):
        """Return the roll angles at the specified time points

        :param time_points: 1D numpy array of length N of absolute time points
        :param interval_indexes: intervals of the time axis where the given time points are expected [optional]

        :return: 1D numpy array of N roll angles
        """
        return self._interpolate_angles(time_points, [Angles.R.value], interval_indexes)

    def __repr__(self):
        axis_str = str(self._time_axis)
        ypr_str = str(self._ypr_angles)
        gso_str = str(self._orbit)

        axis_portion = "Attitude defined on azimuth axis: " + os.linesep + axis_str + os.linesep
        state_vectors_portion = "Yaw Pitch Roll matrix: " + os.linesep + ypr_str + os.linesep
        rotation_order = "Rotation order: {}".format(self.rotation_order.name.upper()) + os.linesep
        reference_frame = "Reference frame: {}".format(self.reference_frame.name.upper()) + os.linesep
        gso_portion = "Attitude info base on orbit:" + os.linesep + gso_str + os.linesep
        return axis_portion + state_vectors_portion + rotation_order + reference_frame + gso_portion


def _check_init_input(gso, azimuth_axis, ypr_matrix, rotation_order, reference_frame):
    _utils.check_type(gso, GeneralSarOrbit, 'gso')
    _utils.check_type(azimuth_axis, are_ax.Axis, 'axis')
    _utils.check_type(ypr_matrix, np.ndarray, 'ypr_matrix')
    _utils.check_ndim_of_numpy_array(ypr_matrix, 2, 'ypr_matrix')
    _utils.check_first_axis_size_of_numpy_array(ypr_matrix, 3, 'ypr_matrix')
    if ypr_matrix.shape[1] < GeneralSarAttitude.get_minimum_number_of_data_points():
        raise RuntimeError(
            "Not enough attitude records provided: {} < {}".format(
                ypr_matrix.shape[1], GeneralSarAttitude.get_minimum_number_of_data_points()))
    if ypr_matrix.shape[1] != azimuth_axis.size:
        raise RuntimeError(
            "Time and attitude records matrix shall have same number of elements: {} != {}".format(
                ypr_matrix.shape[1], azimuth_axis.size))
    _utils.check_dtype_of_numpy_array(ypr_matrix, float, 'ypr_matrix')
    _utils.check_type(rotation_order, str, 'rotation_order')
    _utils.check_type(reference_frame, str, 'reference_frame')


def create_general_sar_attitude(state_vectors: StateVectors, attitude_info: AttitudeInfo) -> GeneralSarAttitude:
    """Create general sar attitude object from state vectors and attitude info metadata

    :param state_vectors: state vectors as a StateVectors metadata object
    :param attitude_info: attitude data as a AttitudeInfo metadata object

    :return: the new GeneralSarAttitude object
    """
    gso = create_general_sar_orbit(state_vectors)
    time_axis = are_ax.RegularAxis((0, attitude_info.time_step, attitude_info.attitude_records_number),
                                   attitude_info.reference_time)
    ypr_matrix = np.vstack((attitude_info.yaw_vector, attitude_info.pitch_vector, attitude_info.roll_vector))
    gsa = GeneralSarAttitude(gso, time_axis, ypr_matrix, attitude_info.rotation_order.value,
                             attitude_info.reference_frame.value)
    return gsa
