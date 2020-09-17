import os
import operator
import logging
import io
import numpy as np
import pyproj as proj
from scipy.interpolate import interp2d
from osgeo import gdal, osr
from arepytools.timing.precisedatetime import PreciseDateTime


def epsg_in_to_epsg_out(xx, yy, zz, epsg_code_in, epsg_code_out):
    p_in = proj.Proj("+init={}".format(epsg_code_in))
    p_out = proj.Proj("+init={}".format(epsg_code_out))
    xx_transf, yy_transf, zz_transf = proj.transform(p_in, p_out, xx, yy, zz)

    return xx_transf, yy_transf, zz_transf


def save_breakpoints(output_folder, breakpoint_names_list, breakpoint_data_list):

    if len(breakpoint_names_list) != len(breakpoint_data_list):
        raise IndexError('Inputs should have same shape')

    for idx, file_name in enumerate(breakpoint_names_list):
        np.save(os.path.join(output_folder, file_name), breakpoint_data_list[idx])


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

            fnf_string = 'TDM_FNF_20_' + lat_tiles[lat_idx] + lon_tiles[lon_idx]
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

    fnf_string_list, geotransform_list, tile_extent_lonlat_list = tandemx_search_fnf_tiles(geographic_boundaries)

    fnf_tile_loaded_list = []
    fnf_loaded_geotransform_list = []

    for tile_idx in np.arange(len(fnf_string_list)):

        fnf_path = os.path.join(fnf_catalogue, fnf_string_list[tile_idx], 'FNF', fnf_string_list[tile_idx] + '.tiff')

        fnf_aux_inf_file_path = os.path.join(
            fnf_catalogue, fnf_string_list[tile_idx], 'AUXFILES', fnf_string_list[tile_idx] + '_INF.txt'
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

                logging.warning('Error: inconsistency for tile' + fnf_string_list[tile_idx] + '\n')

        else:

            logging.warning('Error: tile' + fnf_string_list[tile_idx] + 'not found \n')

        input_image_driver = None

        for idx, line in enumerate(io.open(fnf_aux_inf_file_path, newline='\r\n')):
            month_num = int(line[17:19])
            if month_num == 1:
                month_str = 'JAN'
            elif month_num == 2:
                month_str = 'FEB'
            elif month_num == 3:
                month_str = 'MAR'
            elif month_num == 4:
                month_str = 'APR'
            elif month_num == 5:
                month_str = 'MAY'
            elif month_num == 6:
                month_str = 'JUN'
            elif month_num == 7:
                month_str = 'JUL'
            elif month_num == 8:
                month_str = 'AUG'
            elif month_num == 9:
                month_str = 'SEP'
            elif month_num == 10:
                month_str = 'OCT'
            elif month_num == 11:
                month_str = 'NOV'
            elif month_num == 12:
                month_str = 'DEC'

            date_time = []
            utc_string = line[20:22] + '-' + month_str + '-' + line[12:16] + ' 00:00:00.000000000000'
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

    fnf_string_list, geotransform_list, tile_extent_lonlat_list = tandemx_search_fnf_tiles(
        lon_raster_min, lon_raster_max, lat_raster_min, lat_raster_max
    )

    for tile_idx in np.arange(len(fnf_string_list)):

        fnf_path = os.path.join(out_fnf_path, fnf_string_list[tile_idx], 'FNF', fnf_string_list[tile_idx] + '.tiff')
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
        outdata = driver.Create(fnf_path, Ny, Nx, 1, gdal.GDT_Byte, ['COMPRESS=LZW'])
        outdata.SetGeoTransform(geotransform_list[tile_idx])
        outdata.SetProjection(srs.ExportToWkt())
        outdata.GetRasterBand(1).WriteArray(raster_out)
        outdata.FlushCache()  ##saves to disk!!
        outdata = None
