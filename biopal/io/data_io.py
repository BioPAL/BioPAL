# SPDX-FileCopyrightText: BioPAL <biopal@esa.int>
# SPDX-License-Identifier: MIT

import os
import io
import struct
import logging
import operator
import numpy as np
from scipy.interpolate import interp2d
from matplotlib import pyplot as plt
from biopal.io.xml_io import raster_info
from biopal.utility.constants import EPSG_CODE_LLA
from arepytools.io.productfolder import ProductFolder
from arepytools.timing.precisedatetime import PreciseDateTime
from osgeo import (
    gdal,
    osr,
)
from osgeo.gdalconst import GA_ReadOnly

###############################################################################
# Fields to write and read the Biomass Header in the Binary Files:#############

_STRING_ENCODING = "utf-8"
_date_str_STR_LEN = 33
_date_str_FORMAT = "{}s".format(_date_str_STR_LEN)
_GEO_CORNERS_FORMAT = "ffff"  # lon_min, lon_max, lat_min, lat_max
_UNIQUE_ACQ_ID_STR_LEN = 47
_UNIQUE_ACQ_ID_FORMAT = "{}s".format(_UNIQUE_ACQ_ID_STR_LEN)
_RESOLUTIONS_FORMAT = "ff"  # resolution_m_slant_rg, resolution_m_az
_SENSOR_VELOCITY_FORMAT = "f"  # sensor_velocity
_HEADER_FORMAT = (
    _date_str_FORMAT + _GEO_CORNERS_FORMAT + _UNIQUE_ACQ_ID_FORMAT + _RESOLUTIONS_FORMAT + _SENSOR_VELOCITY_FORMAT
)
_HEADER_FORMAT_SIZE = struct.calcsize(_HEADER_FORMAT)
###############################################################################


def getBiomassHeaderOffsetSize(stack_composition):

    return _HEADER_FORMAT_SIZE


def writeBiomassHeader(
    product_folder,
    channel_idx,
    date_str,
    lon_min,
    lon_max,
    lat_min,
    lat_max,
    unique_acq_id,
    resolution_m_slant_rg,
    resolution_m_az,
    sensor_velocity,
):

    if not isinstance(date_str, str):
        error_message = "date_str should be an Utc string, not a PreciseDateTime object"
        logging.error(error_message)
        raise ValueError(error_message)

    if len(date_str) != _date_str_STR_LEN:
        error_message = "date_str has a different length from an Utc date"
        logging.error(error_message)
        raise ValueError(error_message)

    # encode all the strings with the _STRING_ENCODING format
    encoded_date_str = date_str.encode(_STRING_ENCODING)
    encoded_unique_acq_id = unique_acq_id.encode(_STRING_ENCODING)

    # fill the struct with all data to write
    packed_data = struct.pack(
        _HEADER_FORMAT,
        encoded_date_str,
        lon_min,
        lon_max,
        lat_min,
        lat_max,
        encoded_unique_acq_id,
        resolution_m_slant_rg,
        resolution_m_az,
        sensor_velocity,
    )

    raster_info_read = (
        product_folder.get_channel(channel_idx).metadata.get_metadata_channels(0).get_element("RasterInfo")
    )

    if raster_info_read.header_offset_bytes != _HEADER_FORMAT_SIZE:
        error_message = "Incompatible header offset size, please follow this flow: step 1 -> getBiomassHeaderOffsetSize; step 2 -> append_channel to product_folder with offset from step 1; step 3 -> execute this function"
        logging.error(error_message)
        raise ValueError(error_message)

    raster_file = os.path.join(product_folder.pf_dir_path, raster_info_read.file_name)
    with open(raster_file, "wb") as raster_fid:
        raster_fid.write(packed_data)


def readBiomassHeader(product_folder, channel_idx):

    data_channel_obj = product_folder.get_channel(channel_idx)
    metadata_obj = data_channel_obj.metadata
    metadatachannel_obj = metadata_obj.get_metadata_channels(0)
    ri = metadatachannel_obj.get_element("RasterInfo")
    raster_file = os.path.join(product_folder.pf_dir_path, ri.file_name)

    (
        date_str,
        lon_min,
        lon_max,
        lat_min,
        lat_max,
        unique_acq_id,
        resolution_m_slant_rg,
        resolution_m_az,
        sensor_velocity,
    ) = readBiomassHeader_core(raster_file)

    return (
        date_str,
        lon_min,
        lon_max,
        lat_min,
        lat_max,
        unique_acq_id,
        resolution_m_slant_rg,
        resolution_m_az,
        sensor_velocity,
    )


def readBiomassHeader_core(raster_file):

    # get the product folder RasterInfo and retrive the raster_file name

    # read the raster file (just the header)
    with open(raster_file, "br") as fid:
        packed_data = fid.read(_HEADER_FORMAT_SIZE)

    (
        encoded_date_str,
        lon_min,
        lon_max,
        lat_min,
        lat_max,
        encoded_unique_acq_id,
        resolution_m_slant_rg,
        resolution_m_az,
        sensor_velocity,
    ) = struct.unpack(_HEADER_FORMAT, packed_data)

    date_str = encoded_date_str.decode(_STRING_ENCODING)
    unique_acq_id = encoded_unique_acq_id.decode(_STRING_ENCODING)

    return (
        date_str,
        lon_min,
        lon_max,
        lat_min,
        lat_max,
        unique_acq_id,
        resolution_m_slant_rg,
        resolution_m_az,
        sensor_velocity,
    )


class BiomassL1cRaster:
    """Interfaces wrapper for reading BIOMASS L1c rasters (BioPAL inputs).

    It reads BIOMASS L1c (BioPAL inputs) raster and metadata from path.
    It is suitable to manage any L1c slant range / azumuth BIOMASS data as the BioPAL input SLC stacks, the
    geometric auxiliaries as the ECEFGRID XYZ map, the OffNadirAngles map and so on.

    Parameters
    ----------
    path_dir : str
        path of the directory containing binary data and xml metadata
    channel_to_read : int, default 1
        one-based integer to select the channel to read inside the directory;
        depending on the data, the channel may represent different polarization, baseline, 
        coordinate, etc. as specified in the xml metadata

    Attributes
    ----------
    data : 2D numpy array
        is the loaded data matrix itself
    x_axis : numpy array
        slant range [km] axis
    y_axis : numpy array
        azimuth [km] axis
    x_axis_description : str
        x axis description
    y_axis_description : str
        y axis description

    Returns
    -------
    BiomassL1cRaster : BiomassL1cRaster object

    See Also
    --------
    BiomassL2Raster : used to read biomass L2 and intermediate (geocoded) data
    biopal.utility.plot.plot
    biopal.utility.plot.plot_db
    biopal.utility.plot.plot_abs
    biopal.utility.plot.plot_angle
    biopal.utility.plot.plot_rad2deg : functions to plot BiomassL1cRaster and BiomassL2Raster objects

    Examples
    --------
    >>> data_obj = BiomassL1cRaster(path_dir)
    >>> data_obj = BiomassL1cRaster(path_dir, channel_to_read=1)
    """

    def __init__(self, path_dir, channel_to_read=1):

        self._pf = pf = ProductFolder(path_dir, "r")
        self._data = None  # it is read on the fly, the first time that is requested
        self._channel_to_read = channel_to_read - 1  # zero-based
        # (in the help is described as one-based for consistency with biomassL2raster)
        self._data_type = "SlantRange_Azimuth"  # L1C are the slantrange azimuth slc data

        # axis construction
        data_channel_obj = pf.get_channel(self._channel_to_read)
        metadata_obj = data_channel_obj.metadata
        metadatachannel_obj = metadata_obj.get_metadata_channels(0)
        ri = metadatachannel_obj.get_element("RasterInfo")

        self._x_axis = (ri.samples_start + np.arange(ri.samples) * ri.samples_step) / 1000
        self._y_axis = (np.arange(ri.lines) * ri.lines_step) / 1000
        self._x_axis_description = "slant range [km]"
        self._y_axis_description = "azimuth [km]"

    @property
    def data(self):
        if self._data is None:
            self._data = self._pf.read_data(self._channel_to_read).transpose()
        return self._data

    @property
    def x_axis(self):
        return self._x_axis

    @property
    def y_axis(self):
        return self._y_axis

    @property
    def x_axis_description(self):
        return self._x_axis_description

    @property
    def y_axis_description(self):
        return self._y_axis_description


class BiomassL2Raster:
    """Interfaces wrapper for reading BIOMASS L2 rasters (BioPAL outputs).

    It reads BioPAL L2 geocoded rasrters and metadata from tif file.
    It is suitable to manage:
    BioPAL L2 geocoded east-north output products tif data (default behavior)
    BioPAL intermediate geocoded lat-lon tif data (by setting intermediate_latlon_flag to True)

    Parameters
    ----------
    path_tif : str
        path of the tif file
    band_to_read : int, default 1
        one-based integer to select the tif band to read inside the file;
        depending on the data, multiple bands can be present, as in the AGB estimaton
    intermediate_latlon_flag : bool, default False 
        - False: data axis are considered L2 east/north (default behavior)
        - True:  data axis are considered latitude/longitude (for BioPAL intermediate geocoded products)

    Attributes
    ----------
    data : 2D numpy array
        is the loaded data matrix itself
    x_axis : numpy array
        slant range [km] axis
    y_axis : numpy array
        azimuth [km] axis
    x_axis_description : str
        x axis description
    y_axis_description : str
        y axis description

    Returns
    -------
    BiomassL2Raster : BiomassL2Raster object

    See Also
    --------
    BiomassL1cRaster : used to read biomass L1c data
    biopal.utility.plot.plot
    biopal.utility.plot.plot_db
    biopal.utility.plot.plot_abs
    biopal.utility.plot.plot_angle
    biopal.utility.plot.plot_rad2deg : functions to plot BiomassL1cRaster and BiomassL2Raster objects

    Examples
    --------
    >>> data_obj = BiomassL2Raster(path_tif)

    >>> data_obj = BiomassL2Raster(path_tif, band_to_read=2)

    >>> data_obj = BiomassL2Raster(path_tif, intermediate_latlon_flag=True)

    >>> data_obj = BiomassL2Raster(path_tif, band_to_read=2, intermediate_latlon_flag=True)
    """

    def __init__(self, path_tif, band_to_read=1, intermediate_latlon_flag=None):

        self._data_driver = data_driver = gdal.Open(path_tif, GA_ReadOnly)
        self._data = None  # it is read on the fly, the first time that is requested
        self._band_to_read = band_to_read  # one-based
        self._data_type = "ground"

        # axis construction
        geotransform = data_driver.GetGeoTransform()
        self._x_axis = geotransform[0] + geotransform[1] * np.arange(data_driver.RasterXSize)
        self._y_axis = geotransform[3] + geotransform[5] * np.arange(data_driver.RasterYSize)

        if intermediate_latlon_flag is None:
            # Geocoded EQUI7
            self._x_axis_description = "east [km]"
            self._y_axis_description = "north [km]"
            self._x_axis = self._x_axis / 1000
            self._y_axis = self._y_axis / 1000

        elif intermediate_latlon_flag is True:
            # Geocoded Lat Lon
            self._x_axis_description = "longitude [deg]"
            self._y_axis_description = "latitude [deg]"

    @property
    def data(self):
        if self._data is None:
            self._data = self._data_driver.GetRasterBand(self._band_to_read).ReadAsArray()
        return self._data

    @property
    def x_axis(self):
        return self._x_axis

    @property
    def y_axis(self):
        return self._y_axis

    @property
    def x_axis_description(self):
        return self._x_axis_description

    @property
    def y_axis_description(self):
        return self._y_axis_description


def read_data(folder, pf_name):
    # reads a data:
    # input data is supposed to be an slc in radiometrically calibrated as beta nought or sigma nought
    # it is supposed to contain one or more polarizations (the "SwathInfo" is read to retrive it)
    # it returns a dictionary with keys = polarizations
    # it returns also the dimensions of data

    data_pf_name = os.path.join(folder, pf_name)

    pf = ProductFolder(data_pf_name, "r")
    number_of_pols = pf.get_number_channels()
    data_read = {}
    polid_found = []
    for pol_channel_idx in range(number_of_pols):

        # prepare the metadata elements
        data_channel_obj = pf.get_channel(pol_channel_idx)
        metadata_obj = data_channel_obj.metadata
        metadatachannel_obj = metadata_obj.get_metadata_channels(0)

        # get the ID of the master acquisition:
        di = metadatachannel_obj.get_element("DataSetInfo")
        if not di:
            raise ValueError("data product folder should contain the DataSetInfo to retrive the MASTER ID")
        if di.description.find("Master_swath_") != 0:
            raise ValueError(
                'DataSetInfo description not recognized: it should be a string as "Master_swath_IdOfTheMaster"'
            )
        master_id = di.description[13:]

        # Raster Info
        ri = metadatachannel_obj.get_element("RasterInfo")
        num_samples = ri.samples
        num_lines = ri.lines
        pixel_spacing_slant_rg = ri.samples_step
        pixel_spacing_az = ri.lines_step
        lines_start_utc = str(ri.lines_start)

        # DataSet Info
        di = metadatachannel_obj.get_element("DataSetInfo")
        carrier_frequency_hz = di.fc_hz

        # SwathInfo
        si = metadatachannel_obj.get_element("SwathInfo")
        if not si:
            raise ValueError("data product folder should contain the SwathInfo to retrive the polarization")
        pol_id = si.polarization.name

        polid_found.append(pol_id)
        # Sampling constants
        sc = metadatachannel_obj.get_element("SamplingConstants")
        range_bandwidth_hz = sc.brg_hz

        # hv and vh data (if both are present) are averaged together as (vh+hv)/2"
        if pol_id == "hv" or pol_id == "vh":
            if "vh" in data_read.keys():
                # data (vh or hv) already saved in the dict, add the other data
                data_read["vh"] = (data_read["vh"] + pf.read_data(pol_channel_idx).transpose()) / 2
            else:
                # nor vh nor vv have been saved to dict yet, add first one
                data_read["vh"] = pf.read_data(pol_channel_idx).transpose()
        else:

            data_read[pol_id] = pf.read_data(pol_channel_idx).transpose()

    if len(polid_found) < 3:
        raise ValueError(
            "Input data stack {} should contain #3 or #4 polarizations: hh + hv + vh + vv or hh + hv + vv or hh + vh + vv, found {} instead ".format(
                pf_name, polid_found
            )
        )
    if not "hh" in polid_found or not "vv" in polid_found or not ("hv" in polid_found or "vh" in polid_found):
        raise ValueError(
            "Input data stack {} should contain #3 or #4 polarizations: hh + hv + vh + vv or hh + hv + vv or hh + vh + vv, found {} instead ".format(
                pf_name, polid_found
            )
        )

    return (
        data_read,
        num_samples,
        num_lines,
        pixel_spacing_slant_rg,
        pixel_spacing_az,
        carrier_frequency_hz,
        range_bandwidth_hz,
        master_id,
        lines_start_utc,
    )


def read_auxiliary_single_channel(folder, pf_name):
    # reads Incidence_Angle and Reference_height auxiliary data:
    # it returns a numpy matrix, no dictionaries in this case

    data_pf_name = os.path.join(folder, pf_name)

    if os.path.exists(data_pf_name):
        pf = ProductFolder(data_pf_name, "r")
        number_of_channels = pf.get_number_channels()
        if number_of_channels > 1:
            raise ValueError(
                "Input auxiliary data is supposed to have just one channel, and not # {}".format(number_of_channels)
            )

        aux_read = pf.read_data(0).transpose()
    else:
        aux_read = None
        logging.warning("Path " + data_pf_name + " does not exist.")

    return aux_read


def read_auxiliary_multi_channels(folder, pf_name, valid_acq_id_to_read=None, read_raster_info=False):
    # reads a KZ product:
    # it is supposed to be a pf containing "N" channels, with "N" the number of acquisitions in a stack
    # the acquisition name is read from the SwathInfo "Swath" field
    # it returns a dictionary with keys = acquisition_id ( which is the "Swath")

    data_pf_name = os.path.join(folder, pf_name)

    if os.path.exists(data_pf_name):
        pf = ProductFolder(data_pf_name, "r")
        number_of_acq = pf.get_number_channels()

        data_read = {}
        for channel_idx in range(number_of_acq):

            # prepare the metadata elements
            data_channel_obj = pf.get_channel(channel_idx)
            metadata_obj = data_channel_obj.metadata
            metadatachannel_obj = metadata_obj.get_metadata_channels(0)

            # SwathInfo
            si = metadatachannel_obj.get_element("SwathInfo")
            if not si:
                raise ValueError("Input KZ and off_nadir should contain the SwathInfo to retrive the Swath ID")

            if valid_acq_id_to_read is None or (si.swath in valid_acq_id_to_read):

                data_read[si.swath] = pf.read_data(channel_idx).transpose()

        # Raster Info
        ri = metadatachannel_obj.get_element("RasterInfo")
        num_samples = ri.samples
        num_lines = ri.lines
        pixel_spacing_slant_rg = ri.samples_step
        pixel_spacing_az = ri.lines_step

        raster_info_obj = raster_info(
            num_samples, num_lines, pixel_spacing_slant_rg, pixel_spacing_az, None, None, None, None, None,
        )
    else:
        data_read = None
        raster_info_obj = None
        logging.warning("Path " + data_pf_name + " does not exist.")

    if read_raster_info:
        return data_read, raster_info_obj
    else:
        return data_read


def read_ecef_grid(folder, pf_name):
    # reads an ECEFGRID:
    # it is supposed to be a pf containing exactly 3 channels, with X, Y and Z coordinates
    # the X,Y or Z is got from the DataSetInfo Description field which is supposed
    # to be exactly a string like :Auxiliary data: X ECEF GRID [m] (example for X)
    # it returns a dictionary with keys = coordinate_id (X, Y or Z)

    data_pf_name = os.path.join(folder, pf_name)

    if os.path.exists(data_pf_name):
        pf = ProductFolder(data_pf_name, "r")
        number_of_coords = pf.get_number_channels()
        if not number_of_coords == 3:
            raise ValueError(
                "Input ECEF GRID should contain #3 channels with X,Y and Z coordinates: #{}  channels have been found.".format(
                    number_of_coords
                )
            )

        coordinates_read = {}
        for coord_channel_idx in range(number_of_coords):

            # prepare the metadata elements
            data_channel_obj = pf.get_channel(coord_channel_idx)
            metadata_obj = data_channel_obj.metadata
            metadatachannel_obj = metadata_obj.get_metadata_channels(0)

            # DataSetInfo
            di = metadatachannel_obj.get_element("DataSetInfo")
            if not di:
                raise ValueError("Input ECEF GRID should contain the DataSetInfo to retrive the Description")
            coord_id = di.description[16]
            if not coord_id == "X" and not coord_id == "Y" and not coord_id == "Z":
                raise ValueError(
                    'Cannot retrive coordinate name from DataSetInfo description: description should be a string as: "Auxiliary data: X ECEF GRID [m]", instead it is: "'
                    + di.description
                    + '"'
                )

            coordinates_read[coord_id] = pf.read_data(coord_channel_idx).transpose()

    else:
        coordinates_read = None
        logging.warning("Path " + data_pf_name + " does not exist.")

    return coordinates_read


def tandemx_search_fnf_tiles(geographic_boundaries):
    # geographic_boundaries:
    # is a namedlist with four fields: lon_min, lon_max, lat_min and lat_max

    lon_raster_min = geographic_boundaries.lon_min
    lon_raster_max = geographic_boundaries.lon_max
    lat_raster_min = geographic_boundaries.lat_min
    lat_raster_max = geographic_boundaries.lat_max

    fnf_string_list = []
    geotransform_list = []
    tile_extent_lonlat_list = []

    tile_extent_lat = 1  # deg
    pixel_spacing_lat = 1.8 / 3600  # deg

    lat_start = np.arange(-89, 89, tile_extent_lat) + pixel_spacing_lat / 2

    tile_extent_lon_list = np.zeros(lat_start.shape, dtype=int)
    tile_extent_lon_list[np.logical_and(lat_start >= -89, lat_start < -80)] = 4  # deg
    tile_extent_lon_list[np.logical_and(lat_start >= -80, lat_start < -60)] = 2
    tile_extent_lon_list[np.logical_and(lat_start >= -60, lat_start < 60)] = 1
    tile_extent_lon_list[np.logical_and(lat_start >= 60, lat_start < 80)] = 2
    tile_extent_lon_list[np.logical_and(lat_start >= 80, lat_start < 89)] = 4

    pixel_spacing_lon_list = np.zeros(lat_start.shape, dtype=float)
    pixel_spacing_lon_list[np.logical_and(lat_start >= -89, lat_start < -80)] = 6.4 / 3600  # deg
    pixel_spacing_lon_list[np.logical_and(lat_start >= -80, lat_start < -60)] = 3.6 / 3600
    pixel_spacing_lon_list[np.logical_and(lat_start >= -60, lat_start < 60)] = 1.8 / 3600
    pixel_spacing_lon_list[np.logical_and(lat_start >= 60, lat_start < 80)] = 3.6 / 3600
    pixel_spacing_lon_list[np.logical_and(lat_start >= 80, lat_start < 89)] = 6.4 / 3600

    lat_tiles = ["S{:02d}".format(89 - l) for l in range(0, 89, tile_extent_lat)] + [
        "N{:02d}".format(l) for l in range(0, 89, tile_extent_lat)
    ]

    lat_first_tile = np.max(lat_start[lat_raster_min > lat_start])
    lat_first_tile = np.where(lat_start == lat_first_tile)[0][0]
    lat_last_tile = np.min(lat_start[lat_raster_max <= lat_start])
    lat_last_tile = np.where(lat_start == lat_last_tile)[0][0]

    lat_start = lat_start[lat_first_tile:lat_last_tile]
    tile_extent_lon_list = tile_extent_lon_list[lat_first_tile:lat_last_tile]
    pixel_spacing_lon_list = pixel_spacing_lon_list[lat_first_tile:lat_last_tile]
    lat_tiles = lat_tiles[lat_first_tile:lat_last_tile]

    for lat_idx in np.arange(len(lat_start)):

        pixel_spacing_lon = pixel_spacing_lon_list[lat_idx]
        tile_extent_lon = tile_extent_lon_list[lat_idx]

        lon_start = np.arange(-180, 180, tile_extent_lon) - pixel_spacing_lon / 2

        lon_tiles = ["W{:03d}".format(180 - l) for l in range(0, 180, tile_extent_lon)] + [
            "E{:03d}".format(l) for l in range(0, 180, tile_extent_lon)
        ]

        lon_first_tile = np.max(lon_start[lon_raster_min > lon_start])
        lon_first_tile = np.where(lon_start == lon_first_tile)[0][0]
        lon_last_tile = np.min(lon_start[lon_raster_max <= lon_start])
        lon_last_tile = np.where(lon_start == lon_last_tile)[0][0]

        lon_start = lon_start[lon_first_tile:lon_last_tile]
        lon_tiles = lon_tiles[lon_first_tile:lon_last_tile]

        for lon_idx in np.arange(len(lon_start)):

            fnf_string = "TDM_FNF_20_" + lat_tiles[lat_idx] + lon_tiles[lon_idx]
            geotransform = [
                lon_start[lon_idx],
                pixel_spacing_lon,
                0.0,
                lat_start[lat_idx] + tile_extent_lat,
                0.0,
                -pixel_spacing_lat,
            ]
            tile_extent_lon_lat = [tile_extent_lon, tile_extent_lat]

            fnf_string_list.append(fnf_string)
            geotransform_list.append(geotransform)
            tile_extent_lonlat_list.append(tile_extent_lon_lat)

    return fnf_string_list, geotransform_list, tile_extent_lonlat_list


def tandemx_fnf_read(fnf_catalogue, geographic_boundaries):
    # geographic_boundaries:
    # is a namedlist with four fields: lon_min, lon_max, lat_min and lat_max

    (fnf_string_list, geotransform_list, tile_extent_lonlat_list,) = tandemx_search_fnf_tiles(geographic_boundaries)

    fnf_tile_loaded_list = []
    fnf_loaded_geotransform_list = []

    for tile_idx in np.arange(len(fnf_string_list)):

        fnf_path = os.path.join(fnf_catalogue, fnf_string_list[tile_idx], "FNF", fnf_string_list[tile_idx] + ".tiff",)

        fnf_aux_inf_file_path = os.path.join(
            fnf_catalogue, fnf_string_list[tile_idx], "AUXFILES", fnf_string_list[tile_idx] + "_INF.txt",
        )

        input_image_driver = gdal.Open(fnf_path, 0)

        if input_image_driver is not None:

            Ny = input_image_driver.RasterYSize
            Nx = input_image_driver.RasterXSize
            fnf_mask = input_image_driver.ReadAsArray(0, 0, Nx, Ny)

            fnf_geotransform = input_image_driver.GetGeoTransform()

            map_object = map(operator.sub, list(fnf_geotransform), geotransform_list[tile_idx])
            diff_list = list(map_object)
            values_are_different = [coord_value for coord_value in diff_list if abs(coord_value) > np.finfo(float).eps]
            if not values_are_different:

                fnf_tile_loaded_list.append(fnf_mask)
                fnf_loaded_geotransform_list.append(fnf_geotransform)

            else:

                logging.warning("Error: inconsistency for tile" + fnf_string_list[tile_idx] + "\n")

        else:

            logging.warning("Error: tile" + fnf_string_list[tile_idx] + "not found \n")

        input_image_driver = None

        for idx, line in enumerate(io.open(fnf_aux_inf_file_path, newline="\r\n")):
            month_num = int(line[17:19])
            if month_num == 1:
                month_str = "JAN"
            elif month_num == 2:
                month_str = "FEB"
            elif month_num == 3:
                month_str = "MAR"
            elif month_num == 4:
                month_str = "APR"
            elif month_num == 5:
                month_str = "MAY"
            elif month_num == 6:
                month_str = "JUN"
            elif month_num == 7:
                month_str = "JUL"
            elif month_num == 8:
                month_str = "AUG"
            elif month_num == 9:
                month_str = "SEP"
            elif month_num == 10:
                month_str = "OCT"
            elif month_num == 11:
                month_str = "NOV"
            elif month_num == 12:
                month_str = "DEC"

            date_time = []
            utc_string = line[20:22] + "-" + month_str + "-" + line[12:16] + " 00:00:00.000000000000"
            if idx == 0:
                date_time = PreciseDateTime().set_from_utc_string(utc_string)
            else:
                date_time = max(date_time, PreciseDateTime().set_from_utc_string(utc_string))

    return fnf_tile_loaded_list, date_time, fnf_loaded_geotransform_list


def tandemx_fnf_write(out_fnf_path, fnf_raster, lon_raster, lat_raster):

    lat_raster_min = np.min(lat_raster)
    lon_raster_min = np.min(lon_raster)
    lat_raster_max = np.max(lat_raster)
    lon_raster_max = np.max(lon_raster)

    (fnf_string_list, geotransform_list, tile_extent_lonlat_list,) = tandemx_search_fnf_tiles(
        lon_raster_min, lon_raster_max, lat_raster_min, lat_raster_max
    )

    for tile_idx in np.arange(len(fnf_string_list)):

        fnf_path = os.path.join(out_fnf_path, fnf_string_list[tile_idx], "FNF", fnf_string_list[tile_idx] + ".tiff",)
        directory = os.path.dirname(fnf_path)
        if not os.path.exists(directory):
            os.makedirs(directory)

        lon_out = np.arange(
            geotransform_list[tile_idx][0],
            geotransform_list[tile_idx][0] + tile_extent_lonlat_list[tile_idx][0] + geotransform_list[tile_idx][1],
            geotransform_list[tile_idx][1],
        )
        lat_out = np.arange(
            geotransform_list[tile_idx][3],
            geotransform_list[tile_idx][3] - tile_extent_lonlat_list[tile_idx][1] + geotransform_list[tile_idx][5],
            geotransform_list[tile_idx][5],
        )
        lon_out = lon_out[lon_out <= geotransform_list[tile_idx][0] + tile_extent_lonlat_list[tile_idx][0]]
        lat_out = lat_out[lat_out >= geotransform_list[tile_idx][3] - tile_extent_lonlat_list[tile_idx][1]]

        raster_interp = interp2d(lon_raster, lat_raster, fnf_raster, fill_value=0)
        raster_out = raster_interp(lon_out, lat_out)
        raster_out = np.round(raster_out)

        Nx, Ny = raster_out.shape
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(4326)
        driver = gdal.GetDriverByName("GTiff")
        outdata = driver.Create(fnf_path, Ny, Nx, 1, gdal.GDT_Byte, ["COMPRESS=LZW"])
        outdata.SetGeoTransform(geotransform_list[tile_idx])
        outdata.SetProjection(srs.ExportToWkt())
        outdata.GetRasterBand(1).WriteArray(raster_out)
        outdata.FlushCache()  ##saves to disk!!
        outdata = None


def tiff_formatter(
    data_in, out_fname, geotransform, gdal_data_format, projection=None, multi_layers_tiff=False, time_tag=None,
):

    if ".tiff" in out_fname:
        out_fname = out_fname[0:-5]
    elif ".tif" in out_fname:
        out_fname = out_fname[0:-4]

    if isinstance(data_in, list) and multi_layers_tiff:
        # write multi layer data in same tiff

        if isinstance(geotransform[0], list):
            geotransform = geotransform[0]

        out_tiff_fname = out_fname + ".tif"
        num_layers = len(data_in)
        # formats and saves the input data in GEO-TIFF
        if type(data_in[0]) == str:
            data_temp = np.load(data_in[0])
            Nx, Ny = data_temp.shape
            del data_temp
        else:
            Nx, Ny = data_in[0].shape

        driver = gdal.GetDriverByName("GTiff")

        outdata = driver.Create(out_tiff_fname, Ny, Nx, num_layers, gdal_data_format)
        if time_tag:
            outdata.SetMetadata({"time_tag": time_tag}, "TIFFTAG_DATETIME")

        if projection:
            outdata.SetProjection(projection)

        else:

            srs = osr.SpatialReference()
            srs.ImportFromEPSG(np.int(EPSG_CODE_LLA[5:]))
            outdata.SetProjection(srs.ExportToWkt())

        outdata.SetGeoTransform(geotransform)

        for idx, data in enumerate(data_in):

            if type(data) == str:
                outdata.GetRasterBand(idx + 1).WriteArray(np.load(data))

            else:
                outdata.GetRasterBand(idx + 1).WriteArray(data)

        outdata.FlushCache()  ##saves to disk!!
        outdata = None

    elif isinstance(data_in, list) and not multi_layers_tiff:
        # write each data in a different tiff

        out_tiff_fname = []
        for idx, data in enumerate(data_in):
            out_tiff_fname.append(out_fname + "_fnftile" + str(idx) + ".tif")

            # formats and saves the input data in GEO-TIFF
            Nx, Ny = data.shape

            driver = gdal.GetDriverByName("GTiff")

            outdata = driver.Create(out_tiff_fname[idx], Ny, Nx, 1, gdal_data_format)
            if time_tag:
                outdata.SetMetadata({"time_tag": time_tag}, "TIFFTAG_DATETIME")

            if projection:
                outdata.SetProjection(projection)

            else:
                srs = osr.SpatialReference()
                srs.ImportFromEPSG(np.int(EPSG_CODE_LLA[5:]))
                outdata.SetProjection(srs.ExportToWkt())

            outdata.SetGeoTransform(geotransform[idx])
            outdata.GetRasterBand(1).WriteArray(data)
            outdata.FlushCache()  ##saves to disk!!
            outdata = None

    else:

        # write the single input data to tiff
        out_tiff_fname = out_fname + ".tif"

        # formats and saves the input data in GEO-TIFF
        Nx, Ny = data_in.shape

        driver = gdal.GetDriverByName("GTiff")

        outdata = driver.Create(out_tiff_fname, Ny, Nx, 1, gdal_data_format)
        if time_tag:
            outdata.SetMetadata({"time_tag": time_tag}, "TIFFTAG_DATETIME")

        if projection:
            outdata.SetProjection(projection)

        else:
            srs = osr.SpatialReference()
            srs.ImportFromEPSG(np.int(EPSG_CODE_LLA[5:]))
            outdata.SetProjection(srs.ExportToWkt())

        outdata.SetGeoTransform(geotransform)
        outdata.GetRasterBand(1).WriteArray(data_in)
        outdata.FlushCache()  ##saves to disk!!
        outdata = None

    return out_tiff_fname
