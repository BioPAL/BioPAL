# SPDX-FileCopyrightText: Aresys S.r.l. <info@aresys.it>
# SPDX-License-Identifier: MIT

"""
Constants module
----------------

Example of usage:

.. code-block:: python

    import arepytools.constants as cst
    print(cst.LIGHT_SPEED)
"""

# Speed of light
LIGHT_SPEED = 299792458.0
"""
Speed of light in vacuum (m/s).
"""

# Units of measure
SECOND_STR = 's'
"""
Second symbol.
"""

HERTZ_STR = 'Hz'
"""
Hertz symbol.
"""

JOULE_STR = 'j'
"""
Joule symbol.
"""

DEGREE_STR = 'deg'
"""
Degree symbol.
"""

RAD_STR = 'rad'
"""
Radian symbol.
"""

UTC_STR = 'Utc'
"""
Coordinated Universal Time abbreviation.
"""

METER_STR = 'm'
"""
Meter symbol.
"""

# Metric prefixes
KILO = 1e3
"""
Number of units in one kilounit (kilounits-to-units conversion factor).
"""

MEGA = KILO**2
"""
Number of units in one megaunit (megaunits-to-units conversion factor).
"""

GIGA = KILO**3
"""
Number of units in one gigaunit (gigaunits-to-units conversion factor).
"""

MILLI = KILO**-1
"""
Number of units in one milliunit (milliunits-to-units conversion factor).
"""

MICRO = KILO**-2
"""
Number of units in one microunit (microunits-to-units conversion factor).
"""

NANO = KILO**-3
"""
Number of units in one nanounit (nanounits-to-units conversion factor).
"""

PICO = KILO**-4
"""
Number of units in one picounit (picounits-to-units conversion factor).
"""
