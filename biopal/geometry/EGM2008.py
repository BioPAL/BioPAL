"""
Project: BIODEMPP
Class to retrieve EGM2008 correction

Inputs: 
    db_path (string) (optional): indicates path to geoid database. Requires a local copy of the 
             EGM2008 data-set in GIS format to be downloaded from:
                "http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm2008/egm08_gis.html"
             All uncompressed tiles are to be provided in single directory (db_path)

Outputs: geoid array

Authors: Marc Jaeger & Muriel Pinheiro
DLR. April 2020
"""

# Python libraries
import glob
import os

import numpy as np
import scipy.ndimage as ndi

import gdal
from gdalconst import *


class EGM2008(object):
    def __init__(self, db_path, verbose=False):

        t_files = glob.glob(os.path.join(db_path, "*", "*", "prj.adf"))
        self.tiles = [gdal.Open(f, GA_ReadOnly) for f in t_files]
        self.tform = [t.GetGeoTransform() for t in self.tiles]

        if verbose:
            print("EGM2008: Found %i tiles..." % (len(self.tiles)))

    def get(self, lat, lon):
        geoid = np.zeros(lon.shape, dtype="f4") + np.nan
        for t, f in zip(self.tiles, self.tform):
            lon_px = (lon - f[0]) / f[1]
            lat_px = (lat - f[3]) / f[5]
            mask = (lon_px >= 0) * (lon_px <= t.RasterXSize) * (lat_px >= 0) * (lat_px <= t.RasterYSize)
            vi = np.flatnonzero(mask)
            if vi.size <= 0:
                continue

            tile = np.array(t.GetRasterBand(1).ReadAsArray())
            int_c = np.stack((lat_px.flat[vi], lon_px.flat[vi]))
            geoid.flat[vi] = ndi.map_coordinates(tile, int_c, order=3)

        return geoid
