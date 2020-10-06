import logging
import numpy as np
from scipy.interpolate import griddata
from biopal.utility.utility_functions import epsg_in_to_epsg_out
from biopal.utility.constants import (
    EPSG_CODE_ECEF,
    EPSG_CODE_LLA,
    MAX_EARTH_RADIUS_METERS,
    )


def geocoding_init(ecef_grid, rg_vec_subs, az_vec_subs, min_spacing_m):

    ### Convert ECEFGRID pixel to pixel from X+Y+Z to L+L+H (with for example epsg_in_to_epsg_out): the result will NOT be a uniform grid of lat lon
    logging.info('convert ECEFGRID pixel to pixel from X Y Z to Lat Lon Alt...')
    lon_in, lat_in, alt_in = epsg_in_to_epsg_out(
        ecef_grid['X'][rg_vec_subs, :][:, az_vec_subs],
        ecef_grid['Y'][rg_vec_subs, :][:, az_vec_subs],
        ecef_grid['Z'][rg_vec_subs, :][:, az_vec_subs],
        EPSG_CODE_ECEF,
        EPSG_CODE_LLA,
    )

    lat_lon_step = np.rad2deg(min_spacing_m / MAX_EARTH_RADIUS_METERS)

    logging.info('    geocoding latitude and longitude step used: {} [deg]'.format(lat_lon_step))
    lat_lon_gap = 0  # 1*lat_lon_step

    # latitude should be decrescent
    used_lat_step = -np.abs(lat_lon_step)
    lat_regular_vector = np.arange(np.nanmax(lat_in) - lat_lon_gap, np.nanmin(lat_in) + lat_lon_gap, used_lat_step)
    # longitude should be crescent
    used_lon_step = np.abs(lat_lon_step)
    lon_regular_vector = np.arange(np.nanmin(lon_in) - lat_lon_gap, np.nanmax(lon_in) + lat_lon_gap, used_lon_step)

    lonMeshed_out, latMeshed_out = np.meshgrid(lon_regular_vector, lat_regular_vector)

    valid_values_mask = np.invert(np.isnan(lat_in) + np.isnan(lon_in) + np.isnan(alt_in))

    return (
        lon_in,
        lat_in,
        alt_in,
        lonMeshed_out,
        latMeshed_out,
        valid_values_mask,
        lon_regular_vector,
        used_lon_step,
        lat_regular_vector,
        used_lat_step,
    )


def geocoding(data_to_geocode, lon_in, lat_in, lonMeshed_out, latMeshed_out, valid_values_mask, interp_method='linear'):
    f = griddata(
        np.array([lat_in[valid_values_mask].flatten(), lon_in[valid_values_mask].flatten()]).T,
        data_to_geocode[valid_values_mask].flatten(),
        np.array([latMeshed_out.flatten(), lonMeshed_out.flatten()]).T,
        method=interp_method,
    )  # interpolate to fill holes
    data_ground = f.reshape(latMeshed_out.shape)

    return data_ground