"""
Channel module
------------------------------
"""
import os
from .metadata import MetaData, MetaDataChannel
from ._utils import read_raster, write_raster, write_metadata, data_type_dict_aresys, byte_order_dict_aresys, EOpenMode


class Channel:
    """
    ProductFolder's Channel class
    """

    def __init__(self, file_name, metadata, open_mode):
        self.__open_mode = EOpenMode(open_mode)
        if isinstance(metadata, MetaData):
            self.__metadata = metadata
        else:
            raise TypeError("Unsupported metadata type: {}".format(metadata.__class__.__name__))
        if not (os.path.isfile(file_name)) and self.open_mode == EOpenMode('r'):
            raise OSError(0, "file not found", file_name)
        self.__file_name = file_name

    @classmethod
    def from_raster_info(cls, file_name, raster_info, open_mode):
        metadata = MetaData()
        metadata_channel = MetaDataChannel()
        metadata_channel.insert_element(raster_info)
        metadata.append_channel(metadata_channel)
        return cls(file_name, metadata, open_mode)

    def get_sampling_constants(self, meta_data_ch_index=0):
        """
        Getter of the SamplingConstants

        :param meta_data_ch_index: index of the metadata channel
        :return: :class:`arepytools.io.metadata.SamplingConstants`
        """
        return self.metadata.get_metadata_channels(meta_data_ch_index).get_element('SamplingConstants')

    def get_pulse(self, meta_data_ch_index=0):
        """
        Getter of the Pulse

        :param meta_data_ch_index: index of the metadata channel
        :return: :class:`arepytools.io.metadata.Pulse`
        """
        return self.metadata.get_metadata_channels(meta_data_ch_index).get_element('Pulse')

    def get_raster_info(self, meta_data_ch_index=0):
        """
        Getter of the RasterInfo

        :param meta_data_ch_index: index of the metadata channel
        :return: :class:`arepytools.io.metadata.RasterInfo`
        """
        return self.metadata.get_metadata_channels(meta_data_ch_index).get_element('RasterInfo')

    def get_dataset_info(self, meta_data_ch_index=0):
        """
        Getter of the DataSetInfo

        :param meta_data_ch_index: index of the metadata channel
        :return: :class:`arepytools.io.metadata.DataSetInfo`
        """
        return self.metadata.get_metadata_channels(meta_data_ch_index).get_element('DataSetInfo')

    def get_state_vectors(self, meta_data_ch_index=0):
        """
        Getter of the StateVectors

        :param meta_data_ch_index: index of the metadata channel
        :return: :class:`arepytools.io.metadata.StateVectors`
        """
        return self.metadata.get_metadata_channels(meta_data_ch_index).get_element('StateVectors')

    def get_attitude_info(self, meta_data_ch_index=0):
        """
        Getter of the AttitudeInfo

        :param meta_data_ch_index: index of the metadata channel
        :return: :class:`arepytools.io.metadata.AttitudeInfo`
        """
        return self.metadata.get_metadata_channels(meta_data_ch_index).get_element('AttitudeInfo')

    def get_acquisition_time_line(self, meta_data_ch_index=0):
        """
        Getter of the AcquisitionTimeLine

        :param meta_data_ch_index: index of the metadata channel
        :return: :class:`arepytools.io.metadata.AcquisitionTimeLine`
        """
        return self.metadata.get_metadata_channels(meta_data_ch_index).get_element('AcquisitionTimeLine')

    def get_ground_corner_points(self, meta_data_ch_index=0):
        """
        Getter of the GroundCornerPoints

        :param meta_data_ch_index: index of the metadata channel
        :return: :class:`arepytools.io.metadata.GroundCornerPoints`
        """
        return self.metadata.get_metadata_channels(meta_data_ch_index).get_element('GroundCornerPoints')

    def get_burst_info(self, meta_data_ch_index=0):
        """
        Getter of the BurstInfo

        :param meta_data_ch_index: index of the metadata channel
        :return: :class:`arepytools.io.metadata.BurstInfo`
        """
        return self.metadata.get_metadata_channels(meta_data_ch_index).get_element('BurstInfo')

    def get_doppler_centroid(self, meta_data_ch_index=0):
        """
        Getter of the DopplerCentroidVector

        :param meta_data_ch_index: index of the metadata channel
        :return: :class:`arepytools.io.metadata.DopplerCentroidVector`
        """
        return self.metadata.get_metadata_channels(meta_data_ch_index).get_element('DopplerCentroidVector')

    def get_doppler_rate(self, meta_data_ch_index=0):
        """
        Getter of the DopplerRateVector

        :param meta_data_ch_index: index of the metadata channel
        :return: :class:`arepytools.io.metadata.DopplerRateVector`
        """
        return self.metadata.get_metadata_channels(meta_data_ch_index).get_element('DopplerRateVector')

    def get_tops_azimuth_modulation_rate(self, meta_data_ch_index=0):
        """
        Getter of the TopsAzimuthModulationRateVector

        :param meta_data_ch_index: index of the metadata channel
        :return: :class:`arepytools.io.metadata.TopsAzimuthModulationRateVector`
        """
        return self.metadata.get_metadata_channels(meta_data_ch_index).get_element('TopsAzimuthModulationRateVector')

    def get_slant_to_ground(self, meta_data_ch_index=0):
        """
        Getter of the SlantToGround

        :param meta_data_ch_index: index of the metadata channel
        :return: :class:`arepytools.io.metadata.SlantToGroundVector`
        """
        return self.metadata.get_metadata_channels(meta_data_ch_index).get_element('SlantToGroundVector')

    def get_ground_to_slant(self, meta_data_ch_index=0):
        """
        Getter of the GroundToSlant

        :param meta_data_ch_index: index of the metadata channel
        :return: :class:`arepytools.io.metadata.GroundToSlantVector`
        """
        return self.metadata.get_metadata_channels(meta_data_ch_index).get_element('GroundToSlantVector')

    def get_slant_to_incidence(self, meta_data_ch_index=0):
        """
        Getter of the SlantToIncidence

        :param meta_data_ch_index: index of the metadata channel
        :return: :class:`arepytools.io.metadata.SlantToIncidence`
        """
        return self.metadata.get_metadata_channels(meta_data_ch_index).get_element('SlantToIncidence')

    def get_slant_to_elevation(self, meta_data_ch_index=0):
        """
        Getter of the SlantToElevation

        :param meta_data_ch_index: index of the metadata channel
        :return: :class:`arepytools.io.metadata.SlantToElevation`
        """
        return self.metadata.get_metadata_channels(meta_data_ch_index).get_element('SlantToElevation')

    def get_antenna_info(self, meta_data_ch_index=0):
        """
        Getter of the AntennaInfo

        :param meta_data_ch_index: index of the metadata channel
        :return: :class:`arepytools.io.metadata.AntennaInfo`
        """
        return self.metadata.get_metadata_channels(meta_data_ch_index).get_element('AntennaInfo')

    def get_data_statistics(self, meta_data_ch_index=0):
        """
        Getter of the DataStatistics

        :param meta_data_ch_index: index of the metadata channel
        :return: :class:`arepytools.io.metadata.DataStatistics`
        """
        return self.metadata.get_metadata_channels(meta_data_ch_index).get_element('DataStatistics')

    def get_swath_info(self, meta_data_ch_index=0):
        """
        Getter of the SwathInfo

        :param meta_data_ch_index: index of the metadata channel
        :return: :class:`arepytools.io.metadata.SwathInfo`
        """
        return self.metadata.get_metadata_channels(meta_data_ch_index).get_element('SwathInfo')

    def get_coreg_poly(self, meta_data_ch_index=0):
        """
        Getter of the CoregPolyVector

        :param meta_data_ch_index: index of the metadata channel
        :return: :class:`arepytools.io.metadata.CoregPolyVector`
        """
        return self.metadata.get_metadata_channels(meta_data_ch_index).get_element('CoregPolyVector')

    @property
    def metadata(self):
        return self.__metadata

    @property
    def open_mode(self):
        return self.__open_mode

    @property
    def file_name(self):
        return self.__file_name

    def read_data(self, block_to_read=None):
        """
        Read the raster

        :param block_to_read: portion of the raster to read [lines start, samples start, numbers of lines, number of samples]
        :return: read data
        """

        # Get RasterInfo from metadata
        ri = self.get_raster_info()
        # Get data from RasterInfo
        file_type = data_type_dict_aresys[ri.cell_type.value]
        file_mode = byte_order_dict_aresys[ri.byte_order.value]

        # Read data
        data = read_raster(file_name=self.__file_name, record_len=ri.samples, num_of_records=ri.lines,
                           file_type=file_type, file_mode=file_mode, block_to_read=block_to_read,
                           header_offset=ri.header_offset_bytes, row_prefix=ri.row_prefix_bytes,
                           invalid_value=ri.invalid_value)

        return data

    def read_binary_header(self):
        """
        Read the data header

        :return: data header in bytes
        """
        ri = self.get_raster_info()
        if ri.header_offset_bytes == 0:
            return b''

        with open(self.__file_name, 'rb') as raster_file:
            header = raster_file.read(ri.header_offset_bytes)

        return header

    def write_binary_header(self, header):
        """
        Write the data header

        :param header: header in bytes
        """
        if self.open_mode == EOpenMode('r'):
            raise ReadOnlyChannel

        ri = self.get_raster_info()
        header_size = len(header)
        if ri.header_offset_bytes != header_size:
            raise RuntimeError(
                "Header size incompatible with header offset: {} != {}".format(
                    header_size, ri.header_offset_bytes))

        with open(self.__file_name, 'wb') as raster_file:
            raster_file.write(header)

    def write_data(self, data, start_point=(0, 0)):
        """
        Write the raster

        :param data:
        :param start_point:
        """
        if self.open_mode == EOpenMode('r'):
            raise ReadOnlyChannel

        # Get RasterInfo from metadata
        ri = self.get_raster_info()
        # Get data from RasterInfo
        file_type = data_type_dict_aresys[ri.cell_type.value]
        file_mode = byte_order_dict_aresys[ri.byte_order.value]
        lines_to_write, samples_to_write = data.shape
        max_samples = samples_to_write + start_point[1]
        max_lines = lines_to_write + start_point[0]
        if max_lines <= ri.lines and max_samples <= ri.samples:
            # Write data
            write_raster(file_name=self.__file_name, data=data, record_len=ri.samples, num_of_records=ri.lines,
                         file_type=file_type, file_mode=file_mode, start_point=start_point,
                         header_offset=ri.header_offset_bytes, row_prefix=ri.row_prefix_bytes,
                         invalid_value=ri.invalid_value)
        else:
            raise DataOutOfRasterInfoLimits

    def write_metadata(self):
        """
        Write metadata to xml
        """
        if self.open_mode == EOpenMode('r'):
            raise ReadOnlyChannel
        write_metadata(self.metadata, self.file_name + '.xml')


class ReadOnlyChannel(IOError):
    def __init__(self):
        super().__init__('Operation not allowed: channel open in read mode')


class DataOutOfRasterInfoLimits(IOError):
    def __init__(self):
        super().__init__('The the data to write exceeds the RasterInfo limits')
