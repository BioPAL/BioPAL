"""
Ellipsoid module
----------------
"""

import numpy as np


class Ellipsoid:
    """
    Ellipsoid class
    """

    def __init__(self, first_semi_axis, second_semi_axis):
        """Initialize the object with the specified ellipsoid parameters

        :param first_semi_axis: length of the first semi-axis of the ellipsoid
        :param second_semi_axis: length of the second semi-axis of the ellipsoid
        """

        self._semi_major_axis = float(max(first_semi_axis, second_semi_axis))
        self._semi_minor_axis = float(min(first_semi_axis, second_semi_axis))

        self._semi_axes_ratio_min_max = self._semi_minor_axis / self._semi_major_axis
        self._eccentricity_square = 1 - self._semi_axes_ratio_min_max ** 2
        self._eccentricity = np.sqrt(self._eccentricity_square)
        self._ep2 = 1. / self._semi_axes_ratio_min_max ** 2 - 1

    @property
    def semi_major_axis(self):
        """Semi-major axis length"""
        return self._semi_major_axis

    @property
    def semi_minor_axis(self):
        """Semi-minor axis length"""
        return self._semi_minor_axis

    @property
    def semi_axes_ratio_min_max(self):
        """Min/Max semi-axis ratio"""
        return self._semi_axes_ratio_min_max

    @property
    def eccentricity_square(self):
        """Squared eccentricity"""
        return self._eccentricity_square

    @property
    def eccentricity(self):
        """Eccentricity"""
        return self._eccentricity

    @property
    def ep2(self):
        """Ep2"""
        return self._ep2
