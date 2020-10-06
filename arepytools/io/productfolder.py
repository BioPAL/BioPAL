# SPDX-FileCopyrightText: Aresys S.r.l. <info@aresys.it>
# SPDX-License-Identifier: MIT

"""
ProductFolder module
--------------------
"""
import re
import os
import lxml.etree as etree
from arepytools.io import read_metadata
from ._utils import EOpenMode
from .metadata import RasterInfo
from . import channel

_re_pattern = '_(?!0000)[0-9]{4}'
_header_extension = '.xml'


class ProductFolder:
    """
    ProductFolder class
    """

    def __init__(self, dir_path, open_mode=EOpenMode('r')):

        self.__open_mode = EOpenMode(open_mode)
        self.__pf_dir_path = os.path.abspath(dir_path)
        self.__pf_name = os.path.basename(dir_path)
        self.__channels = list()

        if EOpenMode(open_mode) == EOpenMode('r'):
            if not self.is_productfolder(dir_path):
                raise IsNotAProductFolder
            self.__pf_read(dir_path, open_mode)

        elif EOpenMode(open_mode) == EOpenMode('x') or EOpenMode(open_mode) == EOpenMode('w'):
            if os.path.exists(dir_path):
                raise RuntimeError(
                    "Cannot initialize {} on {}: path already exists".format(self.__class__.__name__, dir_path))

            # Initialize the product folder
            self.__pf_generate(dir_path)
        else:
            raise NotImplementedError("{} mode not supported".format(EOpenMode(open_mode)))

    @property
    def open_mode(self):
        return self.__open_mode

    @property
    def pf_name(self):
        return self.__pf_name

    @property
    def pf_dir_path(self):
        return self.__pf_dir_path

    def __pf_generate(self, dir_path):
        # Create the directory
        os.mkdir(dir_path)
        # Create the manifest
        self.__create_manifest(dir_path)

    @staticmethod
    def __create_manifest(dir_path):
        manifest = etree.Element("AresysProductManifest")
        manifest.set('Version', '2')
        etree.SubElement(manifest, "ProductDescription").text = 'ProductFolder initialized by ArePyTools'
        tree = etree.ElementTree(manifest)
        tree.write(os.path.join(dir_path, 'aresys_product'), pretty_print=True, xml_declaration=True, encoding="utf-8")

    def __pf_read(self, dir_path, open_mode):
        headers = self._get_header_files_in_dir(dir_path)
        raster_extension = self.get_data_file_extension(dir_path)
        for header in headers:
            metadata = read_metadata(os.path.join(dir_path, header))
            raster_filename = header[0:-len(_header_extension)] + raster_extension
            raster_path = os.path.join(dir_path, raster_filename)
            self.__channels.append(channel.Channel(raster_path, metadata, open_mode))

    def append_channel(self, lines, samples, data_type, header_offset=0, row_prefix=0, byte_order='LITTLEENDIAN'):
        """
        Append a channel to the ProductFolder

        :param lines: number of lines
        :param samples: number of samples
        :param data_type: raster data type (see :class:`arepytools.io.metadata.ECellType`)
        :param header_offset: header offset in the raster
        :param row_prefix: row prefix in the raster
        :param byte_order: byte order of the raster (see :class:`arepytools.io.metadata.EByteOrder`)
        """
        if self.open_mode == EOpenMode('r'):
            raise ReadOnlyProductFolder
        # Define the channel file name
        current_num_of_channels = self.get_number_channels()
        _, pf_name = os.path.split(self.pf_dir_path)
        suffix = "_{:0>4}".format(current_num_of_channels + 1)
        chan_file_path = os.path.join(self.pf_dir_path, pf_name + suffix)

        # Create the metadata
        _, file_name = os.path.split(chan_file_path)
        current_raster_info = RasterInfo(lines, samples, data_type, file_name, header_offset, row_prefix, byte_order)

        # Create the new channel
        current_channel = channel.Channel.from_raster_info(chan_file_path, current_raster_info, self.open_mode)

        # Add the channel to the product folder
        self.__channels.append(current_channel)

    @staticmethod
    def _get_header_files_in_dir(dir_path) -> list:
        files = os.listdir(dir_path)
        pf_name = os.path.split(dir_path)[-1]
        header_list = [f for f in files if
                       re.fullmatch("{}".format(pf_name) + _re_pattern + _header_extension, f) is not None]
        header_list.sort()
        return header_list

    @staticmethod
    def _get_raster_files_in_dir(dir_path) -> list:
        files = os.listdir(dir_path)
        pf_name = os.path.split(dir_path)[-1]
        raster_list = [f for f in files if re.fullmatch("{}".format(pf_name) + _re_pattern, f) is not None]
        raster_list.sort()
        return raster_list

    @classmethod
    def is_productfolder(cls, dir_path) -> bool:
        """
        Verify if the folder at the path in input is a ProductFolder

        :param dir_path: path of the folder to check
        :return: True if the folder is a ProductFolder False instead
        """
        if not os.path.isdir(dir_path):
            return False

        files = os.listdir(dir_path)

        # Check manifest existence
        if files.count('aresys_product') != 1:
            return False

        # Exist the raster for each header
        raster_extension = cls.get_data_file_extension(dir_path)
        for header in cls._get_header_files_in_dir(dir_path):
            raster_filename = header[0:-len(_header_extension)] + raster_extension
            if files.count(raster_filename) != 1:
                return False

        # Exist the header for each raster
        for raster_file in cls._get_raster_files_in_dir(dir_path):
            if files.count(raster_file + _header_extension) != 1:
                return False

        return True

    @classmethod
    def get_data_file_extension(cls, dir_path) -> str:
        manifest_filename = os.path.join(dir_path, 'aresys_product')
        manifest_root = etree.parse(manifest_filename).getroot()
        extension_node = manifest_root.find('DataFileExtension')

        if extension_node is not None and extension_node.text is not None:
            return extension_node.text
        return ''

    def get_number_channels(self) -> int:
        """
        Get the number of channels in the ProductFolder

        :return: number of channels
        """
        return len(self.__channels)

    def read_binary_header(self, channel_index):
        """
        Read the header of the raster of a channel of the ProductFolder

        :param channel_index: index of the channel to read
        :return: header in bytes
        """
        return self.__channels[channel_index].read_binary_header()

    def write_binary_header(self, channel_index, header):
        """
        Write data to a channel of the ProductFolder

        :param channel_index: index of the channel to write
        :param header: header in bytes
        """
        if self.open_mode == EOpenMode('r'):
            raise ReadOnlyProductFolder
        self.__channels[channel_index].write_binary_header(header)

    def read_row_prefix(self, channel_index, line_index):
        """
        Read the row prefix of a given line of a channel of the ProductFolder

        :param channel_index: index of the channel
        :param line_index: raster line index
        :return: row prefix as byte sequence
        """
        return self.__channels[channel_index].read_row_prefix(line_index)

    def write_row_prefix(self, channel_index, line_index, row_prefix):
        """
        Write the row prefix of a given line of a channel of the ProductFolder

        :param channel_index: index of the channel
        :param line_index: raster line index
        :param row_prefix: row prefix as bytes sequence
        """
        if self.open_mode == EOpenMode('r'):
            raise ReadOnlyProductFolder
        self.__channels[channel_index].write_row_prefix(line_index, row_prefix)

    def read_data(self, channel_index, block_to_read=None):
        """
        Read the raster of a channel of the ProductFolder

        :param channel_index: index of the channel to read
        :param block_to_read: portion of the raster to read [start line, start sample, number of lines, number of
            samples]
        :return: read data
        """
        data = channel.Channel.read_data(self.__channels[channel_index], block_to_read)
        return data

    def write_data(self, channel_index, data, start_point=(0, 0)):
        """
        Write data to a channel of the ProductFolder

        :param channel_index: index of the channel to write
        :param data: data to write
        :param start_point: coordinates of the first pixel to write
        """
        if self.open_mode == EOpenMode('r'):
            raise ReadOnlyProductFolder
        self.__channels[channel_index].write_data(data, start_point)

    def write_metadata(self, channel_index):
        """
        Write the metadata to the xml for a channel of the ProductFolder

        :param channel_index: index of the channel to write
        """
        if self.open_mode == EOpenMode('r'):
            raise ReadOnlyProductFolder
        self.__channels[channel_index].write_metadata()

    def get_channel(self, channel_index):
        """
        Get a channel of the ProductFolder

        :param channel_index: index of the channel to return
        :return: :class:`arepytools.io.metadata.channel.Channel` -- Channel for the index in input
        """
        return self.__channels[channel_index]


class IsNotAProductFolder(ValueError):
    """Invalid product folder error"""

    def __init__(self):
        super().__init__('The directory in input is not a valid ProductFolder')


class ReadOnlyProductFolder(IOError):
    """Read only product folder error"""

    def __init__(self):
        super().__init__('Operation not allowed: product folder open in read mode')
