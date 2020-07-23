"""
ArePyTools: the Python toolset for SAR data processing
"""

_MIN_PYTHON_VERSION = '3.5'

import sys                                                               # noqa
assert sys.version_info >=\
    tuple((int(v) for v in _MIN_PYTHON_VERSION.split('.'))),\
    "ArePyTools requires Python {} or higher".format(_MIN_PYTHON_VERSION)

import logging                                                           # noqa
logging.getLogger(__name__).addHandler(logging.NullHandler())

__version__ = '1.0.0b2'
