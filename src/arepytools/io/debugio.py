# SPDX-FileCopyrightText: Aresys S.r.l. <info@aresys.it>
# SPDX-License-Identifier: MIT

"""
Binary debug data file management module
----------------------------------------
"""

import struct
import numpy as np

# The debug data file binary header is expected to be a sequence of three 32-bit integers, representing:
# * the cell type identifier
# * the number of samples of the data matrix
# * the number of lines of the data matrix
_HEADER_FORMAT = 'iii'
_HEADER_SIZE = struct.calcsize(_HEADER_FORMAT)

_CELLTYPE_TO_DTYPE_MAP = {
    0: 'float32',
    1: 'complex64',
    2: 'float64',
    3: 'complex128',
    4: 'int16',
    5: 'int8',
    6: 'int32',
    7: 'int64',
}


def read_debug(filename):
    """Read binary debug data file content.

    :param filename: binary debug data file name

    :return: a numpy array containing the data read from the input binary debug data file
    """

    with open(filename, mode='rb') as fdesc:
        header = fdesc.read(_HEADER_SIZE)

        celltype, samples, lines = struct.unpack(_HEADER_FORMAT, header)

        try:
            dtype = _CELLTYPE_TO_DTYPE_MAP[celltype]
        except KeyError as exc:
            raise TypeError('Unknown cell type (cell type id = {})'.format(celltype)) from exc

        if samples == 0 or lines == 0:
            return np.empty((lines, samples), dtype)

        data = np.fromfile(fdesc, dtype)

        return data.reshape((lines, samples))
