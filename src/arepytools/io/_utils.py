# SPDX-FileCopyrightText: Aresys S.r.l. <info@aresys.it>
# SPDX-License-Identifier: MIT

"""
IO utils module
---------------
"""
import os.path
import enum
import numpy as np
import lxml.etree as etree
from . import parsing

data_type_dict = {0: 'f4',
                  1: 'c8',
                  2: 'i2',
                  3: 'i4',
                  4: 'u2',
                  5: 'u4',
                  6: 'b',
                  7: 'B',
                  8: 'i2, i2',
                  9: 'c16',
                  10: 'f8',
                  11: 'i1, i1',
                  }

byte_order_dict = {'ieee-be': '>',
                   'ieee-le': '<',
                   0: '<',
                   }

data_type_dict_aresys = {'FLOAT32': 0,
                         'FLOAT_COMPLEX': 1,
                         'INT16': 2,
                         'INT32': 3,
                         'UINT16': 4,
                         'UINT32': 5,
                         'INT8': 6,
                         'UINT8': 7,
                         'INT16_COMPLEX': 8,
                         'DOUBLE_COMPLEX': 9,
                         'FLOAT64': 10,
                         'INT8_COMPLEX': 11,
                         }

byte_order_dict_aresys = {'BIGENDIAN': 'ieee-be',
                          'LITTLEENDIAN': 'ieee-le',
                          }

_UNSUPPORTED_TYPES = (data_type_dict[8], data_type_dict[11])


class EOpenMode(enum.Enum):
    """Initialization modes

    This enumeration contains all the available initialization modalities
    for the ProductFolder object

    """

    #: open a valid product
    open = 'r'

    #: create a product, if possible (raises an error if it already exists)
    create = 'x'

    #: open a valid product, or, if it not exists, try to create it
    #: (deprecated)
    open_or_create = 'a'

    #: create a product, if possible (raises an error if it already exists)
    create_or_overwrite = 'w'


def read_raster(raster_file_name, num_of_samples, num_of_lines, data_type_id=0,
                binary_ordering_mode=0, block_to_read=None, header_offset=0, row_prefix=0):
    """
    Read a raster file to a numpy array

    :param raster_file_name:
    :param num_of_samples:
    :param num_of_lines:
    :param data_type_id:
    :param binary_ordering_mode:
    :param block_to_read:
    :param header_offset:
    :param row_prefix:
    :return:
    """

    if header_offset < 0:
        raise ValueError("header_offset should be non-negative")

    if row_prefix < 0:
        raise ValueError("row_prefix should be non-negative")

    if data_type_id in data_type_dict.keys():
        data_type = np.dtype(data_type_dict[data_type_id])
    else:
        raise ValueError("Uknown data type id: {}".format(data_type_id))

    if data_type in _UNSUPPORTED_TYPES:
        file_type_string = list(data_type_dict_aresys.keys())[
            list(data_type_dict_aresys.values()).index(data_type_id)]
        raise RuntimeError(
            "Read data from raster of type {} currently not supported.".format(file_type_string))

    file_data_type = np.dtype(
        "{}{}".format(byte_order_dict[binary_ordering_mode], data_type_dict[data_type_id]))

    # Compute the items to read
    if block_to_read is None:
        lines_to_read = num_of_lines
        samples_to_read = num_of_samples
        first_line = 0
        first_sample = 0
    else:
        lines_to_read = block_to_read[2]
        samples_to_read = block_to_read[3]
        first_line = block_to_read[0]
        first_sample = block_to_read[1]

    if first_line < 0:
        raise ValueError("First line to read should be non-negative")

    if first_sample < 0:
        raise ValueError("First sample to read should be non-negative")

    if first_line + lines_to_read > num_of_lines:
        raise ValueError("Block to read exceeds max num lines")

    if first_sample + samples_to_read > num_of_samples:
        raise ValueError("Block to read exceeds max num samples")

    # Read data from file
    with open(raster_file_name, 'rb') as fdesc:
        if samples_to_read == num_of_samples and row_prefix == 0:
            offset_byte = header_offset + first_line * num_of_samples * data_type.itemsize
            data = np.fromfile(fdesc, dtype=file_data_type, count=lines_to_read * samples_to_read,
                               offset=offset_byte)

            return data.reshape((lines_to_read, samples_to_read))

        data = np.empty((lines_to_read, samples_to_read), dtype=data_type)

        offset_byte = (first_line * num_of_samples + first_sample) * data_type.itemsize + \
            row_prefix * (first_line + 1) + header_offset
        fdesc.seek(offset_byte, 0)

        offset_line_byte = (num_of_samples - samples_to_read) * data_type.itemsize + row_prefix

        for line in range(lines_to_read):
            offset_byte = offset_line_byte if line > 0 else 0
            data[line, :] = np.fromfile(fdesc, dtype=file_data_type, count=samples_to_read,
                                        offset=offset_byte)

        return data


def read_metadata(file_name):
    """
    Read the metadata from xml

    :param file_name:
    :return:
    """
    with open(file_name, 'r') as input_file:
        xml_string = input_file.read()

    pyxb_metadata = parsing.pyxb_metadata.CreateFromDocument(xml_string)
    aresys_metadata = parsing.pyxb_aresys_conversion.run_conversion_to_aresys(pyxb_metadata)

    return aresys_metadata


def write_raster(raster_file_name, data, num_of_samples, num_of_lines, data_type_id=0,
                 binary_ordering_mode=0, writing_point=(0, 0), header_offset=0, row_prefix=0):
    """
    Write data to a raster file

    :param raster_file_name:
    :param data:
    :param num_of_samples:
    :param num_of_lines:
    :param data_type_id:
    :param binary_ordering_mode:
    :param writing_point:
    :param header_offset:
    :param row_prefix:
    """

    if header_offset < 0:
        raise ValueError("header_offset should be non-negative")
    if row_prefix < 0:
        raise ValueError("row_prefix should be non-negative")

    if len(writing_point) != 2 or writing_point[0] < 0 or writing_point[1] < 0:
        raise ValueError("Writing point should have two non-negative elements")

    if data_type_dict[data_type_id] in _UNSUPPORTED_TYPES:
        file_type_string = list(data_type_dict_aresys.keys())[
            list(data_type_dict_aresys.values()).index(data_type_id)]
        raise ValueError(
            "Write data to raster of type {} currently not supported.".format(file_type_string))

    # Convert data to data type
    data_type = np.dtype(byte_order_dict[binary_ordering_mode] + data_type_dict[data_type_id])
    data_to_write = np.array(data, dtype=data_type)

    first_line, first_sample = writing_point
    lines_to_write, samples_to_write = data_to_write.shape

    if first_sample + samples_to_write > num_of_samples:
        raise ValueError("Input data exceeds max num samples")

    if first_line + lines_to_write > num_of_lines:
        raise ValueError("Input data exceeds max num lines")

    if os.path.isfile(raster_file_name):
        open_mode = 'r+b'
    else:
        open_mode = 'wb'

    with open(raster_file_name, open_mode) as fdesc:
        raster_size = int(
            header_offset + row_prefix * num_of_lines + data_type.itemsize * num_of_samples * num_of_lines)

        # check raster has the correct size
        fdesc.seek(0, 2)
        file_size = fdesc.tell()
        if file_size < raster_size:
            last_element_position = raster_size - data_type.itemsize
            fdesc.seek(last_element_position)
            fdesc.write(np.array(0, dtype=data_type).tobytes())
        fdesc.seek(0, 2)
        file_size = fdesc.tell()
        if file_size != raster_size:
            raise RuntimeError("Unexpected file size")

        for line in range(0, lines_to_write):
            # Move file cursor to correct position
            absolute_line = line + first_line
            absolute_position = first_sample + num_of_samples * absolute_line
            past_row_prefix_size = row_prefix * (line + first_line + 1)
            write_position = int(
                absolute_position * data_type.itemsize + header_offset + past_row_prefix_size)
            fdesc.seek(write_position)

            # Write
            line_to_write = data_to_write[line].tobytes()
            fdesc.write(line_to_write)

        fdesc.seek(0, 2)
        if fdesc.tell() != raster_size:
            raise RuntimeError("Unexpected final raster size")


def get_line_size(samples, cell_type, row_prefix_size):
    """
    Get the size in bytes of a line in the raster (including row prefix).

    :param samples: number of samples per line
    :param cell_type: data type
    :param row_prefix_size: row prefix size in bytes
    :return: line size in bytes
    """
    data_type_id = data_type_dict_aresys[cell_type.value]
    data_type = np.dtype(data_type_dict[data_type_id])
    return samples * data_type.itemsize + row_prefix_size


def write_metadata(metadata, file_path):
    """
    Write the metadata to a xml

    :param metadata:
    :param file_path:
    """
    metadata_o = parsing.pyxb_aresys_conversion.run_conversion_pyxb(metadata)
    xml_string = metadata_o.toxml(encoding="utf-8", element_name='AresysXmlDoc')

    root = etree.fromstring(xml_string)
    tree = etree.ElementTree(root)
    tree.write(file_path, pretty_print=True, xml_declaration=True, encoding="utf-8")
