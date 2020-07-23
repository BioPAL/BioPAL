"""
Coordinate conversions module
-----------------------------
"""

from typing import Union
import numpy as np
from arepytools.geometry.wgs84 import WGS84
from arepytools import _utils

_NUMBER_OF_ITERATIONS = 5

_A_MAX_SQUARE_MINUS_A_MIN_SQUARE_DIVIDED_BY_A_MAX = \
    (WGS84.semi_major_axis ** 2 - WGS84.semi_minor_axis ** 2) / WGS84.semi_major_axis
_A_MIN_TIMES_EP2 = WGS84.semi_minor_axis * WGS84.ep2


def llh2xyz(coordinates: Union[list, np.ndarray]) -> np.ndarray:
    """
    Conversion from llh coordinates to xyz.

    :param coordinates: a numpy ndarray of shape (3,), (3, 1) or (3, N) or a list.
    :return: a two dimensional numpy array of shape (3, N) with xyz coordinates
    """
    coordinates = _utils.input_data_to_numpy_array_with_checks(coordinates, dtype=float, first_axis_size=3,
                                                               name="coordinates")
    if coordinates.ndim == 1:
        coordinates = coordinates.copy()
        coordinates.shape = (3, 1)

    lat = coordinates[0, :]
    lon = coordinates[1, :]
    h = coordinates[2, :]

    big_n = WGS84.semi_major_axis / np.sqrt(1 - WGS84.eccentricity_square * np.sin(lat) ** 2)
    r = (big_n + h) * np.cos(lat)

    x = r * np.cos(lon)
    y = r * np.sin(lon)
    z = (big_n * WGS84.semi_axes_ratio_min_max ** 2 + h) * np.sin(lat)

    return np.vstack((x, y, z)).astype(float)


def xyz2llh(coordinates: Union[list, np.ndarray]) -> np.ndarray:
    """
    Conversion from xyz coordinates to llh.

    :param coordinates: a numpy ndarray of shape (3,), (3, 1) or (3, N) or a list.
    :return: a two dimensional numpy array of shape (3, N) with llh coordinates
    """
    coordinates = _utils.input_data_to_numpy_array_with_checks(coordinates, dtype=float, first_axis_size=3,
                                                               name="coordinates")
    if coordinates.ndim == 1:
        coordinates = coordinates.copy()
        coordinates.shape = (3, 1)

    x = coordinates[0, :]
    y = coordinates[1, :]
    z = coordinates[2, :]

    lon = np.arctan2(y, x)

    k = np.maximum(np.absolute(x), np.absolute(y))
    k[np.where(k == 0)] = 1.

    c1 = x / k
    c2 = y / k

    r = k * np.sqrt(c1 ** 2 + c2 ** 2)

    beta = np.arctan2(z, r * WGS84.semi_axes_ratio_min_max)
    sin_beta, cos_beta = np.sin(beta), np.cos(beta)

    lat = np.arctan2(z + _A_MIN_TIMES_EP2 * sin_beta ** 3,
                     r - _A_MAX_SQUARE_MINUS_A_MIN_SQUARE_DIVIDED_BY_A_MAX * cos_beta ** 3)
    beta_new = np.arctan2(WGS84.semi_axes_ratio_min_max * np.sin(lat), np.cos(lat))

    for _ in range(_NUMBER_OF_ITERATIONS):
        if np.array_equal(beta, beta_new):
            break

        beta = beta_new
        sin_beta, cos_beta = np.sin(beta), np.cos(beta)

        lat = np.arctan2(z + _A_MIN_TIMES_EP2 * sin_beta ** 3,
                         r - _A_MAX_SQUARE_MINUS_A_MIN_SQUARE_DIVIDED_BY_A_MAX * cos_beta ** 3)
        beta_new = np.arctan2(WGS84.semi_axes_ratio_min_max * np.sin(lat), np.cos(lat))

    # Compute height
    sin_phi = np.sin(lat)
    big_n = WGS84.semi_major_axis / np.sqrt(1 - sin_phi ** 2 * WGS84.eccentricity_square)
    h = r * np.cos(lat) + (z + WGS84.eccentricity_square * big_n * sin_phi) * sin_phi - big_n

    return np.vstack((lat, lon, h)).astype(float)
