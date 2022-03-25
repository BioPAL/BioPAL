# SPDX-FileCopyrightText: BioPAL <biopal@esa.int>
# SPDX-License-Identifier: MIT

import numpy as np
from scipy.interpolate import interp1d
from biopal.io.xml_io import raster_info


def apply_calibration_screens(data_slc, slc_raster_info, cal_screens, screens_raster_info, master_id):

    for idx, pf_name in enumerate(data_slc.keys()):

        if idx == 0:
            rg_axis_index = 0
            az_axis_index = 1

            # input original axis
            max_val_in_rg = screens_raster_info.num_samples * screens_raster_info.pixel_spacing_slant_rg
            rg_ax_in = np.arange(0, max_val_in_rg, raster_info.pixel_spacing_slant_rg)

            max_val_in_az = screens_raster_info.num_lines * screens_raster_info.pixel_spacing_az
            az_ax_in = np.arange(0, max_val_in_az, screens_raster_info.pixel_spacing_az)

            # output interpolated axis
            max_val_out_rg = slc_raster_info.num_samples * slc_raster_info.pixel_spacing_slant_rg
            rg_ax_out = np.arange(0, max_val_out_rg, slc_raster_info.pixel_spacing_slant_rg)

            max_val_out_az = slc_raster_info.num_lines * slc_raster_info.pixel_spacing_az
            az_ax_out = np.arange(0, max_val_out_az, slc_raster_info.pixel_spacing_az)

        # interpolate current screen over the data definiton axes:
        cal_scree_interp = cal_screens[pf_name]

        # interpolation of data along axis_index
        interp_fun = interp1d(rg_ax_in, cal_scree_interp, rg_axis_index, bounds_error=0)
        cal_scree_interp = interp_fun(rg_ax_out)
        interp_fun = interp1d(az_ax_in, cal_scree_interp, az_axis_index, bounds_error=0)
        cal_scree_interp = interp_fun(az_ax_out)

        for pol_id in data_slc[pf_name]:
            data_slc[pf_name][pol_id] = np.multiply(data_slc[pf_name][pol_id], np.exp(-1j * cal_scree_interp))

    return data_slc
