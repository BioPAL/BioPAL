Getting Started
===============

Section overview:

* requirements
* installation procedure
* processor setup and run procedure

for a complete guide in setting-up and runnning BioPAL, also have a look at the :ref:`Tutorials` section.

Requirements
------------

Python >= 3.7.1

The needed python packages are specified in the file `requirements`_.

.. _requirements: https://github.com/BioPAL/BioPAL/blob/main/requirements.txt


System requirements
-------------------
CPU, no restrictions; suggested at least an Intel Core i5 10th gen 64bit processor or equivanent.

20 GB of free RAM

Windows 10 or Linux, 64bit


Installation
------------

Installation procedure described here makes use of the open-source package management system `conda`_.

.. _conda: https://docs.conda.io/projects/conda/en/latest/

Note that the installation and processor run procedures are different in case of developers and basic users. 
The differences are underlined when needed. 


Installation prerequisites
^^^^^^^^^^^^^^^^^^^^^^^^^^

* `conda <https://docs.conda.io/projects/conda/en/latest/>`_  should be already installed
* In case you are a developer
    * `git <https://git-scm.com/>`_ , should be already installed
    * `tortoisegit <https://tortoisegit.org>`_ , a git GUI for Windows, optional installation  
    * python IDE (i.e. `spyder <https://www.spyder-ide.org/>`_ ,  `vs code <https://code.visualstudio.com>`_),  optional installation


BioPAL installation default option: "pip install" (users)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

BioPAL will be automatically downoladed from PyPI.

Open a command window with *conda* available and follow this procedure.

Create an empty biopal environment: you can customize the *biopal* environment name::

    conda create --name biopal python==3.7.1

Install GDAL library::

    conda activate biopal
    conda install -c conda-forge GDAL=3.5

Install the package::	

    pip install biopal


BioPAL installation for developers only: "pip install -e"
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The code will be editable, thanks to the "-e" option.

First, you need to fork your own copy of BioPAL from web interface at `github BioPAL <https://github.com/BioPAL/BioPAL>`_.
Your private fork url will be something like *https://github.com/your_name_here/BioPAL*

Than open a command window with *conda* available and follow this procedure.

Make a local clone inside an empty folder (*your_installation_folder/*)::
    
    cd your_installation_folder/
    git clone --branch <branchname> <remote-repo-url> .

* *remote-repo-url* is the url created during the fork from `github BioPAL <https://github.com/BioPAL/BioPAL>`_
* *branchname* is the branch to be cloned: if is a new fork there is only a branch called *main*

Create an empty biopal environment: you can customize the *biopal* environment name::

    conda create --name biopal python==3.7.1
		
Install GDAL library::

    conda activate biopal
    conda install -c conda-forge GDAL=3.5

Install the package from BioPAL/ folder::

    cd your_installation_folder/BioPAL/
    pip install -e .


BioPAL datasets
^^^^^^^^^^^^^^^

BioPAL gives easy access to several datasets that are used for examples in the documentation and testing. 
These datasets are hosted on our FTP server and must be downloaded for use. 

Contact <biopal@esa.int> to receive access to the dataset and for more information.


Setup Configuration
^^^^^^^^^^^^^^^^^^^

Quick Start
"""""""""""
Quick start, after installation, is required to get editable and usable biopal xml Input and Configuration files.

Open a command window with *conda* available and follow this procedure.

Quick Start command::

    biopal-quickstart FOLDER

"FOLDER" is the path where usable and editable versions of `Input_File.xml` 
and `Configuration_File.xml` files will be generated.

Run the processor
^^^^^^^^^^^^^^^^^

From FOLDER (see *quickstart* section), open the *Input_File.xml* and verify/update following sections:

* *output_specification->output_folder*: output folder, each run corresponds to a sub-folder formatted with the current date time
* *dataset_query->L1C_repository*: path of the *dataSet* folder with the stacks to be processed
* *dataset_query->auxiliary_products_folder*: path of the *auxiliary_data_pf* folder with parameters related to the data stacks specified in *L1C_repository*
* *dataset_query->L1C_date* and *->geographic_boundaries_polygon* : those fields are already filled with default values ready to be used with the currently available demo dataSets from ESA.

IMPORTANT: all the paths in the Input_File.xml should be ABSOLUTE paths

NOTE: Sample data (L1C_repository dataSets) and auxiliaries (auxiliary_products_folder) can be obtained by writing to <biopal@esa.int>.

Set *Configuration_File.xml* present in FOLDER (see *quickstart*), as desired:
the AGB, FH, FD, TOMO_FH configuration sections have ready default configuration parameters.

Open a command window with *conda* available and follow this procedure.

Activate the biopal environment::
    
    conda activate biopal

Run BioPAL::

    biopal --conf conf_folder inputfilexml

* *inputfilexml*: path of the BioPAL xml input file 
* *conf_folder*:  path of the folder containing BioPAL xml configuration file

*Input_File.xml* and *conf_folder* may be the ones present in FOLDER (generated during *quickstart*), 
or any other custom ones.

Or Run BioPAL with default configurations::

    biopal inputfilexml

Default configurations are equal to the ones generated during *quickstart*.

Or show BioPAL help::

    biopal -h


Run the processor for developers, with a script for debug
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""

How to run the processor with a script to be launced from an IDE.

Create a new *.py* script file as::

    from pathlib import Path
    import sys
    import os
    biopal_path = Path( 'your_installation_folder/BioPAL' )
    sys.path.append( str(biopal_path) )
    os.chdir(biopal_path)
    from biopal.__main__ import biomassL2_processor_run
    input_file_xml_path = biopal_path.joinpath('Input_File.xml')
    conf_folder = 'yourConFolder/'
    biomassL2_processor_run(input_file_xml_path, conf_folder )

*your_installation_folder/BioPAL* is the folder where BioPAL has been git-cloned.
*Input_File.xml* and *conf_folder* may be the ones generated during *quickstart*, or any other custom ones.

Execute the script within your preferred IDE options (i.e. run, debug, breakpoints enabled...).

Read the :doc:`tutorials` section for other scripts, and for manual execution of a BioPAL chain, step by step.

	
GDAL paths troubleshooting
""""""""""""""""""""""""""
The BioPAL GDAL paths are automatically found by the processor after a correct installation procedure.

In case of problems or for particular user cases, it is possible to manually specify such paths
by editing the *Configuration_File.xml* (from FOLDER):

uncomment the *gdal* section and insert your absolute paths for

* *gdal_path*: this is the folder containing the GDAL executables, usually in the */bin* subfolder of GDAL environment (containing e.g., *gdalwarp*, *gdal_translate*,... )
* *gdal_enviroment_path*: this is the GDAL_DATA environment variable path

IMPORTANT: all the gdal paths, if specified in the Configuration_File.xml, should be ABSOLUTE paths

TIP: the above paths depend on your machine environment. 

GDAL has been automatically installed during the above procedure of conda environment creation; 
for a standard installation with conda, the paths should be found in paths similar to the following (where *xxx* is an alphanumeric string depending on the GDAL version installed)

Windows:

* gdal_path (i.e.): *C:\ProgramData\Anaconda3\pkgs\libgdal-xxx\Library\bin*
* gdal_enviroment_path (i.e.): *C:\ProgramData\Anaconda3\pkgs\libgdal-xxx\Library\share\gdal*

Linux:

* gdal_path (i.e.): */home/user/.conda/envs/biopal/bin*
* gdal_enviroment_path (i.e.): */home/user/.conda/pkgs/libgdal-xxx/share/gdal*