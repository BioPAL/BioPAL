from . import parsing
import os.path
import numpy as np
import enum
import lxml.etree as etree

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

    #: try to create a product, if possible, or, if it exists, try to
    #: overwrite it
    create_or_overwrite = 'w'


def read_raster(file_name, record_len, num_of_records, file_type=0, file_mode=0, block_to_read=None, header_offset=0,
                row_prefix=0, invalid_value=None):
    """
    Read a raster file to a numpy array

    :param file_name:
    :param record_len:
    :param num_of_records:
    :param file_type:
    :param file_mode:
    :param block_to_read:
    :param header_offset:
    :param row_prefix:
    :param invalid_value:
    :return:
    """

    # Check the byte order
    byte_order_str = byte_order_dict[file_mode]

    # Check the data type
    if file_type in data_type_dict.keys():
        data_type = np.dtype(data_type_dict[file_type])
    else:
        raise ValueError

    if data_type in _UNSUPPORTED_TYPES:
        file_type_string = list(data_type_dict_aresys.keys())[list(data_type_dict_aresys.values()).index(file_type)]
        raise RuntimeError(
            "Read data from raster of type {} currently not supported.".format(file_type_string))

    # Compute the items to read
    if block_to_read is None:
        lines_to_read = num_of_records
        samples_to_read = record_len
        first_line = 0
        first_sample = 0
    else:
        lines_to_read = block_to_read[2]
        samples_to_read = block_to_read[3]
        first_line = block_to_read[0]
        first_sample = block_to_read[1]

    # Read data from file
    data = np.empty(shape=(lines_to_read, samples_to_read), dtype=data_type)
    with open(file_name, 'rb') as fd:
        current_pos_byte = 0
        for line in range(0, lines_to_read):
            offset_byte = ((line + first_line) * record_len + first_sample) * data_type.itemsize \
                          - current_pos_byte + row_prefix * (line + first_line + 1) + header_offset
            data[line] = np.core.records.fromfile(fd, formats=data_type, shape=samples_to_read, offset=offset_byte,
                                                  byteorder=byte_order_str)
            current_pos_byte += offset_byte + samples_to_read * data_type.itemsize

    return data


def read_metadata(file_name):
    """
    Read the metadata from xml

    :param file_name:
    :return:
    """
    if os.path.isfile(file_name):
        with open(file_name, 'r') as f:
            xml_string = f.read()
    else:
        raise OSError(0, "file not found", file_name)

    pyxb_metadata = parsing.pyxb_metadata.CreateFromDocument(xml_string)
    aresys_metadata = parsing.pyxb_aresys_conversion.run_conversion_to_aresys(pyxb_metadata)

    return aresys_metadata


def write_raster(file_name, data, record_len, num_of_records=0, file_type=0, file_mode=0, start_point=None,
                 header_offset=0, row_prefix=0, invalid_value=None):
    """
    Write data to a raster file

    :param file_name:
    :param data:
    :param record_len:
    :param num_of_records:
    :param file_type:
    :param file_mode:
    :param start_point:
    :param header_offset:
    :param row_prefix:
    :param invalid_value:
    """

    if data_type_dict[file_type] in _UNSUPPORTED_TYPES:
        file_type_string = list(data_type_dict_aresys.keys())[list(data_type_dict_aresys.values()).index(file_type)]
        raise RuntimeError(
            "Write data to raster of type {} currently not supported.".format(file_type_string))

    # Convert data to data type
    data_type = np.dtype(byte_order_dict[file_mode] + data_type_dict[file_type])
    data_to_write = np.array(data, dtype=data_type)

    # Compute the items to read
    lines_to_write, samples_to_write = data_to_write.shape
    if start_point is None:
        first_line = 0
        first_sample = 0
        record_len = max(samples_to_write, record_len)
    else:
        first_line = start_point[0]
        first_sample = start_point[1]
        record_len = max(first_sample + samples_to_write, record_len)

    max_lines = max(num_of_records, lines_to_write + first_line)
    # Write data to raster
    if os.path.isfile(file_name):
        open_mode = 'r+b'
    else:
        open_mode = 'wb'

    with open(file_name, open_mode) as fd:
        offset_byte = int((record_len * max_lines - 1) * data_type.itemsize + header_offset + row_prefix * max_lines)
        fd.seek(offset_byte)
        fd.write(np.array(0, dtype=data_type).tobytes())

        for line in range(0, lines_to_write):
            # Convert data to bytes
            line_to_write = data_to_write[line].tobytes()
            # Compute offset
            offset_byte = int((record_len *
                               (line + first_line) + first_sample)
                              * data_type.itemsize + header_offset + row_prefix * (line + first_line + 1))

            fd.seek(offset_byte)
            # Write
            fd.write(line_to_write)

        fd.seek(0, 2)
        num_of_records = int(fd.tell() / record_len / data_type.itemsize)
    return record_len, num_of_records


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
