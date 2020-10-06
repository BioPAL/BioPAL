# SPDX-FileCopyrightText: Aresys S.r.l. <info@aresys.it>
# SPDX-License-Identifier: MIT

"""
WGS84 module
------------
"""

from arepytools.geometry._ellipsoid import Ellipsoid

_A_MAX = 6.378137e6  # semi-major axis
_A_MIN = 6.356752314245e6  # semi-minor axis

WGS84 = Ellipsoid(_A_MAX, _A_MIN)
