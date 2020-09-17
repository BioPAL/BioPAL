import os
import sys

###SETTING ENVIROMENT

# append l2 processor code to python path
src_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(src_path)

# set current directory to "l2processor" folder, which contains the main module "main_biomassL2_processor.py"
os.chdir(os.path.join(src_path))

# INSTALLATION_FOLDER is the folder containing "conf", "doc", "input", "src", and"output" subfolders
INSTALLATION_FOLDER = os.path.dirname(src_path)
from biomassL2.main_biomassL2_processor import biomassL2_processor_main

###

### USER settings and Inputs:
# setting gdal paths, this is an example from Anaconda enviroment
gdal_path = r'gdal/bin'
gdal_enviroment_path = r'pkgs/libgdal-2.3.3-h2e7e64b_0/share/gdal'

# Input file:
input_file_xml = os.path.join(INSTALLATION_FOLDER, 'conf/Input_File.xml')
###

# run processor:
biomassL2_processor_main(input_file_xml, gdal_path, gdal_enviroment_path, INSTALLATION_FOLDER)
