"""
Geometry interpolator module
----------------------------
"""

import numpy as np
from arepytools import _utils


def check_eval_function_input(f):
    """
    Interpolator eval function input checker
    """

    def decorated_fun(self, *args_f, **kwargs_f):
        interpolation_axis = _utils.input_data_to_numpy_array_with_checks(args_f[0], ndim=1,
                                                                          name='Interpolation axis')
        return f(self, interpolation_axis, *args_f[1:], **kwargs_f)

    return decorated_fun


class GeometryInterpolator:
    """
    Internal class for handling interpolation in GeneralSarOrbit and in GeneralSarAttitude
    It interpolates vector-valued functions defined on a given axis
    """

    _POLYNOMIAL_ORDER = 4
    _DATA_POINTS_PER_POLY = _POLYNOMIAL_ORDER + 1
    _POLY_CENTER_REL_INTERVAL_INDEX = _DATA_POINTS_PER_POLY // 2

    @classmethod
    def get_min_num_of_data_points(cls):
        """Return the minimum number of data points needed by the interpolator

        :return: minimum number of data points"""
        return cls._DATA_POINTS_PER_POLY

    @classmethod
    def get_polynomial_order(cls):
        """Return the order of the interpolation polynomial

        :return: order of the interpolation polynomial
        """
        return cls._POLYNOMIAL_ORDER

    def __init__(self, data_axis, data):
        """Initialize the object as interpolator of the input data defined on the specified axis

        :param data_axis: the axis of length N on which the data is defined
        :param data: data defined as MxN array of M-dimensional values (e.g. state vectors or attitude angles)
        """
        self._data_axis = data_axis
        self._num_of_data_components = data.shape[0]
        self._num_of_polynomials = self._data_axis.size - self._POLYNOMIAL_ORDER

        # Matrix of polynomials
        self._polynomials_matrix = np.ndarray((self._num_of_polynomials, self._num_of_data_components), dtype=np.poly1d)

        # Relative azimuth axis [s]
        relative_axis = self._data_axis.get_relative_array()
        for poly_index in range(self._num_of_polynomials):
            # Each polynomial domain is defined by the poly_axis
            poly_axis = relative_axis[poly_index: poly_index + self._DATA_POINTS_PER_POLY]

            # Each polynomial relative axis starts at the start of the interval:
            # P(t) = c0 + c1 * (t-poly_time_axis[0]) + ...
            vandermonde_matrix = np.vander((poly_axis - poly_axis[0]).astype(float))

            # Computation of polynomial coefficients
            for component_index in range(self._num_of_data_components):
                values = data[component_index, poly_index: poly_index + self._DATA_POINTS_PER_POLY]
                coefficients = np.linalg.solve(vandermonde_matrix, values).tolist()
                self._polynomials_matrix[poly_index, component_index] = np.poly1d(coefficients)

    @check_eval_function_input
    def eval(self, interpolation_axis, interval_indexes, components_to_interpolate=None):
        """Interpolate the data on the given interpolation axis

        :param interpolation_axis: interpolation axis
        :param interval_indexes: vector or single scalar or None. Denotes the intervals where the given interpolation
        points are expected
        :param components_to_interpolate: list of the components to interpolate, identified by their
        index [optional]

        :return: CxP array of interpolated values, where C is the number of interpolated components and P is the length
        of the interpolation axis
        """
        return self._eval(interpolation_axis, interval_indexes, components_to_interpolate)

    @check_eval_function_input
    def eval_first_derivative(self, interpolation_axis, interval_indexes, components_to_interpolate=None):
        """Interpolate the first derivative of the data on the given interpolation axis

        :param interpolation_axis: interpolation axis
        :param interval_indexes: vector or single scalar or None. Denotes the interval where the given interpolation
        points are expected
        :param components_to_interpolate: list of the components to interpolate, identified by their index [optional]

        :return: CxP array of interpolated values, where C is the number of interpolated components and P is the length
        of the interpolation axis
        """
        return self._eval(interpolation_axis, interval_indexes, components_to_interpolate, np.polyder)

    @check_eval_function_input
    def eval_second_derivative(self, interpolation_axis, interval_indexes, components_to_interpolate=None):
        """Interpolate the second derivative of the data on the given interpolation axis

        :param interpolation_axis: interpolation axis
        :param interval_indexes: vector or single scalar or None. Denotes the interval where the given interpolation
        points are expected
        :param components_to_interpolate: list of the components to interpolate, identified by their index [optional]

        :return: CxP array of interpolated values, where C is the number of interpolated components and P is the length
        of the interpolation axis
        """

        def second_der(p):
            return np.polyder(p, 2)

        return self._eval(interpolation_axis, interval_indexes, components_to_interpolate, second_der)

    def _eval(self, interpolation_axis, interval_indexes, components_to_interpolate=None,
              polynomial_modifier=lambda x: x):
        """Interpolate the data on the given interpolation axis by evaluating interpolation
        polynomials (or a their modified version) on the interpolation axis

        :param interpolation_axis: interpolation axis
        :param interval_indexes: vector or single scalar or None. Denotes the interval where the given time is expected
        :param components_to_interpolate: list of the components to interpolate, identified by their index [optional]
        :param polynomial_modifier: modifier to apply to the interpolation polynomials before evaluation

        :return: CxP array of interpolated values, where C is the number of interpolated components and P is the length
        of the interpolation axis
        """

        if components_to_interpolate is None:
            components_to_interpolate = range(self._num_of_data_components)

        poly_indexes = self.get_polynomial_index(interpolation_points=interpolation_axis,
                                                 interval_indexes=interval_indexes)

        interpolation_results = np.zeros((len(components_to_interpolate), interpolation_axis.size))
        data_axis_array = self._data_axis.get_array()

        for interpolation_point_index, (poly_index, interpolation_point) in enumerate(
                zip(poly_indexes, interpolation_axis)):
            interpolation_point_poly_rel = interpolation_point - data_axis_array[poly_index]

            for out_component_index, component_index in enumerate(components_to_interpolate):
                interpolation_results[out_component_index, interpolation_point_index] = \
                    polynomial_modifier(self._polynomials_matrix[poly_index, component_index])(
                        interpolation_point_poly_rel)

        return interpolation_results

    def get_interval_index(self, interpolation_points):
        """Get the indexes of the intervals where the interpolation_points lie

        :param interpolation_points: vector of required interpolation points

        :return: array of interval indexes
        """
        return self._data_axis.get_interval_id(interpolation_points)

    def get_polynomial_index(self, interpolation_points=None, interval_indexes=None):
        """Get the indexes of the polynomials to be used for interpolation

        :param interpolation_points: vector of required interpolation points
        :param interval_indexes: interval indexes of the interpolation points

        :return: array of polynomial indexes

        If interpolation points are specified, the positions on the interpolation axis are computed internally to find
        the corresponding polynomial indexes;
        If the interval indexes of the interpolation points on the axis are given, the polynomial indexes are computed
        directly (ignoring the interpolation points).
        """
        if interval_indexes is None:
            if interpolation_points is not None:
                interval_indexes = self._data_axis.get_interval_id(interpolation_points)
            else:
                raise RuntimeError("Specify either interpolation points or interval indexes")
        if isinstance(interval_indexes, int):
            interval_indexes = np.asarray([interval_indexes])
        return np.clip(interval_indexes - self._POLY_CENTER_REL_INTERVAL_INDEX, 0, self._num_of_polynomials - 1)
