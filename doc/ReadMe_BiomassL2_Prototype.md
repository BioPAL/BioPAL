# BiomassL2 Prototype Processor, AGB only, 24 JULY 2020

# How to install

The software requirements are:
Operating System: Linux x64

## The package does not need installation, it should be executed under a python enviroment as here described.

Python 3.7 environment, with following packages:

* Equi7Grid 0.0.10
* Pytileproj 0.0.12 ( this is a sub-dependence of the above Equi7Grid )
* GeoPandas 0.7.0
* GDAL 2.3.3
* lxmx 4.4.1
* Matplotlib 2.2.2
* Namedlist 1.7
* Numpy 1.15.4
* Numpydoc 0.8.0
* Packaging 19.0
* Pandas 1.0.5
* Pyproj 2.1.1
* PyXB 1.2.6
* Scipy 1.4.1
* Shapely 1.6.4.post2


## Usage, Preparation of the package:

NOTE: 
the "input" folder with data is needed to perform a run, and is provided separately.

Open the conf/Input_File.xml:
Insert your current directory locations with absolute paths for:
* L1cRepository: this should point over the input/dataSet folder ("input" folder here not present)
* AuxiliaryProductsFolder: this should point over the input/auxiliary_data_pf  folder ("input" folder here not present)
* OutputFolder: is the folder where the output will be generated (note that at each launch a sub-folder with the current date is generated, so any folder will never be overwritten, as better explained in section 5)

Open the src/run_biomassL2_processor_AGB.py

Insert, in the opened py file, your current directory locations with absolute paths for:
* input_file_xml: is the path of the conf/Input_File.xml file
* gdal_path: this is the folder containing the GDAL executables, usually in the "/bin" sub-folder of GDAL environment ( containing e.g. "gdalwarp", "gdal_translate", ... )
* gdal_enviroment_path: this is the GDAL_DATA environment variable path

## Usage, run the processor:
To execute the procesor, run the following script in your python environment:
src/run_biomassL2_processor_AGB.py

For more detalis on processing configuration, refer to the user manual:
ARE-017082_BIOMASS_L2_User_Manual.pdf,  version 1.1
