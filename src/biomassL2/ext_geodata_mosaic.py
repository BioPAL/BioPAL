'''
Project: BIODEMPP
Class to retrieve a DEM or FNF from a local database

Inputs: 
    latmin,latmax,lonmin,lonmax (float): indicagte min and max latitudes/longitudes, in degrees; 
    database_dir (string): indicates path to database; 
    data_type (string): indicates type of data being read 
          supported datatypes are: Copernicus_DEM, TDM_FNF, TDM90_DEM
    output_dir (string): indicates path to output folder which contains Aresys product folder
    geoid_dir (string) (optional): indicates path to geoid database. Requires a local copy of the 
             EGM2008 data-set in GIS format to be downloaded from:
                "http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm2008/egm08_gis.html"
             All uncompressed tiles are to be provided in single directory (geoid_dir)


Expected folder structure of supported data-bases (in principle, it is consistent with tree obtained after unziping):
    CopernicusDEM (example path to database: '/BIODEMPP/CopernicusDEM/'):
    /BIODEMPP/CopernicusDEM/DEM1_SAR_DGE_90_20101225T170840_20131115T170933_ADS_000000_2798.DEM$ ls
        Copernicus_DSM_30_N48_00_E010_00  GSC#CR#ESA#COP-DEM_GLO-90-DGED#20191202#061113.xml  INSPIRE.xml
    TanDEM-X DEM (example path to database: '/BIODEMPP/TDM-DEM/''):    
    /BIODEMPP/TDM-DEM/TDM1_DEM__30_N57W075_V02_C$ ls
        DEM  demProduct.xsd  generalHeader.xsd  TDM1_DEM__30_N57W075.xml  types_inc.xsd
    TanDEM-X FNF (example path to database: '/BIODEMPP/TDM-FNF/'):     
    /BIODEMPP/TDM-FNF/TDM_FNF_20_N41E065$ ls
        FNF  TDM_FNF_20_N41E065_QL.png       

Outputs: mosaicked/cropped data written in Aresys format. Product folder given by
        output_dir+data_type+'_MOS'

Authors: Nestor Yague-Martinez, Jose Luis Bueso Bello, Joel Amao & Muriel Pinheiro
DLR. February 2020
'''

# Python libraries
import os
import shutil
import getpass
import logging
import gdal
import csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET
from xml.dom import minidom

# Internal libraries
# from src.miscellanea.EGM2008 import EGM2008
from biomassL2.EGM2008 import EGM2008
from arepytools.io import metadata
from arepytools.io.productfolder import ProductFolder
from arepytools.io.productfolder import EOpenMode

# from src.inout.XmlIO import add_param
from biomassL2.XmlIO import add_param


class ext_geodata_mosaic:
    def __init__(self, latmin, latmax, lonmin, lonmax, database_dir, data_type, output_dir, geoid_dir=None):

        # log=logging.getLogger('ext_geodata_mosaic')

        # User inputs
        self.input_dir = database_dir
        self.output_dir = output_dir
        self.src = data_type
        self.latmin = latmin
        self.lonmin = lonmin
        self.latmax = latmax
        self.lonmax = lonmax
        self.geoid_dir = geoid_dir

        # create name of Aresys Product folder, Xml and auxiliary XML
        proudct_tag = data_type + '_MOS'
        self.geodata_output_dir = output_dir + proudct_tag
        self.geodata_output_raster = self.geodata_output_dir + '/' + proudct_tag + '_0001'
        self.geodata_output_Xml = self.geodata_output_dir + '/' + proudct_tag + '_0001.xml'
        self.geodata_output_auxXml = self.geodata_output_dir + '/' + proudct_tag + '_auxXml.xml'

        # set invalid value according to data type
        if self.src == 'TDM_FNF':
            self.dtype_string = 'int'
            self.dtype = metadata.ECellType.int8
            self.invalid_value = 0
        else:
            self.dtype_string = 'float'
            self.dtype = metadata.ECellType.float32
            self.invalid_value = -32767.0

        # The CSV Database file is generated only once (if file not available)
        self.DB_output_filename = 'DEMDB_coord.csv'

        # Generate CSV output for database, if non-existent
        self.generate_csv_with_coordinates()

    def __del__(self):
        del self.out_arr

    # -------------------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------------------
    # generate_csv_with_coordinates
    # -------------------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------------------

    def generate_csv_with_coordinates(self):
        log = logging.getLogger('ext_geodata_mosaic')

        self.outputFile_coord = os.path.join(self.input_dir, self.DB_output_filename)
        if not os.path.exists(self.outputFile_coord):
            log.info('Extracting tif information for reference database')
            # files = os.listdir(input_dir)
            file1 = []
            for root, dirs, files in os.walk(self.input_dir):
                for file in files:
                    if self.src == 'Copernicus_DSM':
                        if ('DEM.tif' in file) and (self.src in file):
                            log.info('dem file found %s' % file)
                            file1.append(os.path.join(root, file))
                    elif self.src == 'TDM_FNF':
                        if ('.tif' in file) and (self.src in file):
                            log.info('FNF file found %s' % file)
                            file1.append(os.path.join(root, file))
                    elif self.src == 'TDM90_DEM':
                        if ('DEM.tif' in file) and (self.src in file):
                            log.info('dem file found %s' % file)
                            file1.append(os.path.join(root, file))
                    else:
                        raise ValueError(
                            'Invalid data base name. Valid databases are: Copernicus_DSM, TDM_FNF, TDM_DEM'
                        )

            files = sorted(set(file1))

            output_data = []
            header = ['MINLAT', 'MINLON', 'MAXLAT', 'MAXLON', 'LINK']
            output_data.append(header)
            for file1, i in zip(files, range(len(files))):
                if (self.src in file1) and (r'.tif' in file1):  # Filename: fromglc10v01_0_10.tif
                    log.info('Reading file %s (%i from %i)' % (file1, i, len(files)))
                    file1 = os.path.join(self.input_dir, file1)
                    # Get coordinates
                    gdal1 = gdal.Open(file1)
                    maxLat = gdal1.GetGeoTransform()[3]  # data.geoPoint[0]
                    minLon = gdal1.GetGeoTransform()[0]  # data.geoPoint[1]
                    samLat = -gdal1.GetGeoTransform()[5]  # -data.sampling[0]
                    samLon = gdal1.GetGeoTransform()[1]  # data.sampling[1]
                    minLat = maxLat - (gdal1.RasterYSize) * samLat
                    maxLon = minLon + (gdal1.RasterXSize) * samLon
                    del gdal1

                    output_data.append([minLat, minLon, maxLat, maxLon, file1])

            outputFile = os.path.join(self.input_dir, self.DB_output_filename)
            log.info('Generating output file: %s' % outputFile)
            with open(outputFile, 'w', newline='') as f:
                c = csv.writer(f, delimiter=';')
                c.writerows(output_data)
        else:
            log.info('reference database already available')

        return

    # -------------------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------------------
    # correct for EGM2008
    #
    # -------------------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------------------
    def corr_geoid(self, lat, lon):

        log = logging.getLogger('ext_geodata_mosaic')
        log.info('Correcting geoid')

        try:
            latm = np.tile(lat, (len(lon), 1)).transpose()
            lonm = np.tile(lon, (len(lat), 1))
            geoid = EGM2008(self.geoid_dir, verbose=False).get(latm, lonm)
            geoid = np.nan_to_num(geoid, nan=0)
            # adding geoid height to have heights over the ellipsoid
            self.out_arr += geoid
        except:
            log.warning('Not possible to load geoid! Result corresponds to geoidal heights!')

        return

    # -------------------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------------------
    # get_mosaic
    # -------------------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------------------
    def get_mosaic(self):

        log = logging.getLogger('ext_geodata_mosaic')

        log.info(
            'Obtaining geodata defined by geographical corners: ('
            + str(self.latmin)
            + ','
            + str(self.lonmin)
            + ') and ('
            + str(self.latmax)
            + ','
            + str(self.lonmax)
            + ') deg'
        )

        # Get matching tiles with input coordinates
        df = pd.read_csv(self.outputFile_coord, sep=';')

        minlats = df['MINLAT']
        maxlats = df['MAXLAT']
        minlons = df['MINLON']
        maxlons = df['MAXLON']

        case1 = minlats > self.latmax
        case2 = maxlats < self.latmin
        case3 = minlons > self.lonmax
        case4 = maxlons < self.lonmin
        case = case1 | case2 | case3 | case4  # case1 + case2 + case3 + case4
        indices = np.where(case == 0)[0]  # indices are NOT refered to original dataframe!!!
        links = df['LINK']
        ref_files = np.array(links)[indices].tolist()
        dst_filename = os.path.join(self.input_dir, 'out_dem.tif')

        skip_gdalTrans = True
        if skip_gdalTrans and (len(ref_files) > 0):
            tmp_ds = gdal.Open(ref_files[0])
            xRes = tmp_ds.GetGeoTransform()[1]
            yRes = -tmp_ds.GetGeoTransform()[5]
            out_ds = gdal.Warp(
                dst_filename,
                ref_files,
                format='GTiff',
                outputBounds=(self.lonmin, self.latmin - yRes, self.lonmax + xRes, self.latmax),
            )
        else:
            tmp_ds = gdal.Warp(
                'temp', ref_files, format='MEM'
            )  # gdalWarp needed to mosaic input files. Gdal.translate takes only a single input file as input
            out_ds = gdal.Translate(dst_filename, tmp_ds, projWin=[self.lonmin, self.latmax, self.lonmax, self.latmin])
            # width = lonsize, height = latsize) # Output size = defined input size

        self.out_arr = out_ds.ReadAsArray()  # output DEM
        # rcleanup objects
        del tmp_ds, out_ds

        return

    # -------------------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------------------
    # save_mosaic
    # -------------------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------------------
    def save_mosaic(self, plot=False, forcePermission=True):
        log = logging.getLogger('ext_geodata_mosaic')

        # get geo information
        gdalObj = gdal.Open(self.input_dir + 'out_dem.tif')
        GT = gdalObj.GetGeoTransform()
        nLat, nLon = self.out_arr.shape
        lat = GT[3] + np.arange(nLat) * GT[5]
        lon = GT[0] + np.arange(nLon) * GT[1]
        minLon = lon[0]
        minLat = lat[-1]
        samLat = -GT[5]
        samLon = GT[1]
        del gdalObj

        # correct geoid in case of Copernicus DEM
        if self.src == 'Copernicus_DSM':
            if self.geoid_dir:
                self.corr_geoid(lat, lon)
            else:
                log.warning('Path to EGM2008 not given! Result corresponds to geoidal heights!')

        # remove tempory file
        if os.path.exists(self.input_dir + 'out_dem.tif'):
            os.remove(self.input_dir + 'out_dem.tif')

        # 1. Save raster to disk according to aresys format
        log.info('Saving data and metadata to disk: ' + self.geodata_output_dir)

        # Create output folder
        open_mode = EOpenMode.create_or_overwrite
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
        # Cleanup ProductFolder to avoid crash
        if os.path.exists(self.geodata_output_dir):
            shutil.rmtree(self.geodata_output_dir)
        # Create the ProductFolder
        pf = ProductFolder(self.geodata_output_dir, open_mode)

        # Write to disk
        data_channel_index = 0
        if self.src == 'TDM_FNF':
            self.out_arr = np.asarray(self.out_arr, dtype='int8')
        else:
            self.out_arr = np.asarray(self.out_arr, dtype='float32')
        pf.append_channel(self.out_arr.shape[0], self.out_arr.shape[1], self.dtype)
        pf.write_data(data_channel_index, self.out_arr, (0, 0))

        # 2. Write the Aresys metadata to disk
        chan = pf.get_channel(data_channel_index)
        md = chan.metadata
        md_chan = md.get_metadata_channels(data_channel_index)
        ri = md_chan.get_element('RasterInfo')
        # mapping latitude in lines and longitude in samples
        ri.set_lines_axis(minLat, 'deg', samLat, 'deg')
        ri.set_samples_axis(minLon, 'deg', samLon, 'deg')
        pf.write_metadata(data_channel_index)

        if plot:
            plt.imshow(self.out_arr)
            plt.show()

        # 3. Auxiliary Biodempp metadata (at the moment only storing database, data_path and invalid)
        log.info('Generating BIODEMPP metadata...')
        self.save_biodempp_metadata()

        # 4. Changing permission
        if forcePermission:
            os.system('chmod 664 ' + self.geodata_output_auxXml)
            os.system('chmod 664 ' + self.geodata_output_Xml)
            os.system('chmod 664 ' + self.geodata_output_raster)

        return

    # -------------------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------------------
    # save auxiliary biodempp metadata in auxiliary Xml file (not necessary at the moment)
    # -------------------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------------------
    def save_biodempp_metadata(self):

        log = logging.getLogger('ext_geodata_mosaic')
        log.info('Creating biodempp metadata')

        # Creating an XML file from scratch
        root = ET.Element("Geodata_BIODEMPP")
        geodata = ET.Element("object", name="Auxiliary Parameters")

        # Type of data
        type_data = add_param(
            geodata,
            name="geodata_type",
            unit_text="string",
            datatype_text="string",
            remark_text="",
            value_text=self.src,
        )
        # Invalid value
        inval_value = add_param(
            geodata,
            name="invalid_value",
            unit_text="1",
            datatype_text=self.dtype_string,
            remark_text="",
            value_text=str(self.invalid_value),
        )

        root.append(geodata)
        dom = minidom.parseString(ET.tostring(root))

        # Writing to file
        with open(self.geodata_output_auxXml, 'wb') as outfile:
            outfile.write(dom.toprettyxml(encoding="UTF-8"))

        return

    # -------------------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------------------
    #  run all required steps
    # -------------------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------------------

    def Run(self):

        # ---- Configuration logging
        logging.basicConfig(
            filename='logging.log',
            level='INFO',
            format='%(asctime)s - %(levelname)s - %(name)s: %(message)s',
            datefmt='%d/%m/%Y %I:%M:%S %p',
            filemode='w',
        )
        # log also to console
        stream_handler = logging.StreamHandler()
        # add formatter to stream handler
        stream_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(name)s: %(message)s'))
        logging.getLogger().addHandler(stream_handler)
        # create logger for main script
        log = logging.getLogger('main')
        # --------

        log.info('Starting external geodata mosaicker by user: ' + getpass.getuser())

        self.get_mosaic()

        self.save_mosaic()

        log.info('External geodata mosaicker succesfully finished')

        return
