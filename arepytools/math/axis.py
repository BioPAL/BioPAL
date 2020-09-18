# SPDX-FileCopyrightText: Aresys S.r.l. <info@aresys.it>
# SPDX-License-Identifier: MIT

"""
Axis module
-----------
"""

import os
import typing
import numpy as np
import scipy.interpolate

from arepytools import _utils
from arepytools.timing.precisedatetime import PreciseDateTime

# Types
RealNumber = typing.Union[int, float]
AxisStartType = typing.Union[RealNumber, PreciseDateTime]
RegularAxisTuple = typing.Tuple[RealNumber, RealNumber, int]
GeneralAxisVector = np.ndarray


class Axis:
    """Axis class"""

    def __init__(self, relative_axis: GeneralAxisVector, origin: AxisStartType = 0):
        """
        Constructs an axis in the form origin + relative_axis that starts at origin + relative_axis[0]

        Relative axis should be a monotone array either increasing or decreasing

        :param relative_axis: a numpy array of real numbers
        :param origin: either a real number or a  PreciseDataTime [optional, default 0]
        """
        valid, error_message = _validate_general_axis(relative_axis)
        if not valid:
            raise ValueError(error_message)

        self._relative_axis = relative_axis

        if not isinstance(origin, (PreciseDateTime, int, float)):
            raise ValueError('origin should be either a PreciseDateTime or a real number')

        self._origin = origin

        self._increasing = _is_increasing(self._relative_axis)

    @property
    def start(self) -> AxisStartType:
        """
        :return: origin + relative_axis[0]

        see also origin
        """
        return self._origin + self._relative_axis[0]

    @property
    def origin(self) -> AxisStartType:
        """
        :return: origin

        see also start
        """
        return self._origin

    @property
    def increasing(self) -> bool:
        """
        :return: whether the axis is monotone increasing
        """
        return self._increasing

    @property
    def decreasing(self) -> bool:
        """
        :return: whether the axis is monotone decreasing
        """
        return not self.increasing

    @property
    def mean_step(self) -> RealNumber:
        """
        :return: the average step
        """
        return np.diff(self._relative_axis).mean()

    @property
    def size(self) -> int:
        """Number of axis points"""
        return self._relative_axis.size

    @property
    def length(self) -> RealNumber:
        """
        :return length of the internal part of the axis (N - 1) * step
        """
        return abs(self._relative_axis[-1] - self._relative_axis[0])

    def get_array(self, start: int = 0, stop: int = None) -> np.ndarray:
        """
        Get the underlying axis (default) or a portion of it ([start, stop)) in the form of a numpy.ndarray.

        :param start: start index (included)
        :param stop: stop index (excluded)
        :return: the array
        """
        return self.get_relative_array(start, stop) + self._origin

    def get_relative_array(self, start=0, stop=None) -> np.ndarray:
        """
        Get the underlying relative axis (default) or a portion of it ([start, stop)) in the form of a numpy.ndarray.

        The result is relative to the origin, i.e. the axis = relative_array + origin.

        :param start: start index (included)
        :param stop: stop index (excluded)
        :return: an array of floats
        """
        stop = self.size if stop is None else stop
        _check_range(start, stop, self.size)
        return self._relative_axis[start: stop]

    def get_interval_id(self, values) -> np.ndarray:
        """
        Find the interval of the axis in which each values lies.

        If values lie outside the axis the closest intervals are returned.

        :param values: values to locate on the axis
        :return: np.ndarray with the position of the values in the axis
        """
        return _get_interval_id_not_regular_real_axis(self._relative_axis, self.increasing, values - self._origin)

    def get_interval_id_from_relative(self, values) -> np.ndarray:
        """
        Find the interval of the axis in which each values lies.
        Here values are relative to origin.

        :param values: values to locate on the axis
        :return: np.ndarray with the position of the values in the axis
        """
        return self.get_interval_id(values + self._origin)

    def interpolate(self, position):
        """
        For a given position it returns the value of the axis via linear interpolation.

        It coincides with accessing the array at the given position when the position is an integer,
        otherwise it interpolates the closest values

        :param position: the position (integer of float)
        :return: the corresponding axis point
        """
        return scipy.interpolate.interp1d(range(self.size), self._relative_axis)(position) + self._origin

    def __repr__(self):
        """
        :return: string representation of the object
        """
        axis_repr = "{classname} -- direction: {direction}:" + \
                    "{newline} start: {start}, end: {end}" + \
                    "{newline} size: {size}, length: {length}"
        return axis_repr.format(start=self.start, size=self.size,
                                end=self.start + self.length if self.increasing else -self.length,
                                newline=os.linesep + '--', classname=self.__class__.__name__,
                                direction="increasing" if self.increasing else "decreasing", length=self.length)


class RegularAxis(Axis):
    """Regular axis class"""

    def __init__(self, relative_axis: RegularAxisTuple, origin: AxisStartType = 0):
        """
        A RegularAxis is an Axis where the relative_axis is given in the form:
            (start [float, int], step[float, int], size [int)

        Step can be negative

        :param relative_axis: a tuple (start [float, int], step[float, int], size [int)
        :param origin: either a real number or a PreciseDataTime [optional, default 0]
        """
        # input validation
        valid, error_message = _validate_uniform_axis(relative_axis)
        if not valid:
            raise ValueError(error_message)

        # relative axis initialization
        start, step, size = relative_axis
        super(RegularAxis, self).__init__(np.asarray([start + step * k for k in range(size)]), origin)

    @property
    def step(self) -> RealNumber:
        """Axis step"""
        return self.mean_step

    def interpolate(self, position):
        """
        For a given position it returns the value of the axis via linear interpolation.

        out = position * step + start

        :param position: the position (integer of float)
        :return: the corresponding axis point
        """
        return position * self.step + self.start

    def get_interval_id(self, values) -> np.ndarray:
        """
        Find the interval of the axis in which each values lies.

        If values lie outside the axis the closest intervals are returned.

        :param values: values to locate on the axis
        :return: np.ndarray with the position of the values in the axis
        """
        return _get_interval_id_regular_real_axis(self._relative_axis[0], self.step, self.size, values - self._origin)


def _check_range(start, stop, size):
    if stop > size:
        raise ValueError("Stop ({stop}) should be smaller equal than size: ({size})".format(stop=stop, size=size))
    if start < 0:
        raise ValueError("Start ({start}) should be positive".format(start=start))


def _is_increasing(axis):
    steps = np.diff(axis)
    if (steps > 0).all():
        return True
    elif (steps < 0).all():
        return False
    else:
        raise RuntimeError("Expecting monotone axis")


def _get_interval_id_regular_real_axis(start: RealNumber, step: RealNumber, size: int, values) -> np.ndarray:
    """
    The interval id is intended as the interval starting at the returned position, containing the point.
    If the point is exactly an edge of the interval, the interval starting at the point should be returned.

    :param start: the start of the axis
    :param step: the step of the axis
    :param size: the size of the axis
    :param values: the values to be located in the array
    :return: np.ndarray of the same size as values
    """

    def to_int_and_clip(value):
        return np.clip(int(np.floor(value)), 0, size - 1)

    val = (values - start) / step

    if isinstance(val, np.ndarray):
        return np.array([to_int_and_clip(v) for v in val])
    else:
        return np.array([to_int_and_clip(val)])


def _get_interval_id_not_regular_real_axis(array: np.ndarray, increasing: bool, values: np.ndarray) -> np.ndarray:
    """
    The interval id is intended as the interval starting at the returned position, containing the point.
    If the point is exactly an edge of the interval, the interval starting at the point should be returned.

    :param array
    :param values: the values to be located in the array
    :return: np.ndarray of the same size as values
    """

    values = _utils.input_data_to_numpy_array_with_checks(values, ndim=1)
    out = np.empty((values.size,))
    for k, v in enumerate(values):
        closest_pos = np.argmin(np.abs((array - v)))
        if increasing:
            if array[closest_pos] > v:
                out[k] = closest_pos - 1
            else:
                out[k] = closest_pos
        else:
            if array[closest_pos] < v:
                out[k] = closest_pos - 1
            else:
                out[k] = closest_pos
    return np.clip(out.astype(int), 0, array.size - 1)


def _validate_uniform_axis(uniform_axis: tuple) -> typing.Tuple[bool, str]:
    if not isinstance(uniform_axis, tuple):
        return False, "relative axis should be a tuple"

    if not len(uniform_axis) == 3:  # start, step and size
        return False, "relative axis should be a tuple of size three (start, step, size)"

    if not isinstance(uniform_axis[0], (float, int)):  # start
        return False, "start of relative axis {} should be a real number".format(type(uniform_axis[0]))

    if not isinstance(uniform_axis[1], (float, int)):  # step
        return False, "Step of relative axis {} should be a real number".format(type(uniform_axis[1]))

    if not isinstance(uniform_axis[2], int):  # step
        return False, "Type of relative axis {} should be an integer".format(type(uniform_axis[2]))

    if uniform_axis[2] < 1:  # step
        return False, "Size of relative axis {} should be an integer".format(uniform_axis[2])

    return True, ""


def _validate_general_axis(general_axis: np.ndarray) -> typing.Tuple[bool, str]:
    if not isinstance(general_axis, np.ndarray):
        return False, "relative axis should be a numpy.ndarray"

    if not ((general_axis.dtype == int) or (general_axis.dtype == float)):
        return False, "relative axis dtype {} should be a numpy.ndarray of real numbers".format(general_axis.dtype)

    if not general_axis.size == general_axis.shape[0]:
        return False, 'relative axis should be a vector not a matrix {} != {}'.format(general_axis.size,
                                                                                      general_axis.shape[0])

    return True, ""
