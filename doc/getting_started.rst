Getting Started
===============

In this section requirements and installation procedure of BioPAL are described;
It is also briefly described how to run the processor:
for a complete guide in executing BioPAL and producing results have a look at the Tutorial section.

Requirements
------------

Python 3.7.1 is a minimum requirement. 
The packages required are specified in the file `requirements`_.

.. _requirements: https://github.com/BioPAL/BioPAL/blob/main/requirements.txt


Installation
------------

Installation procedure described here makes use of the open-source package management system `conda`_.

.. _conda: https://docs.conda.io/projects/conda/en/latest/

Note that the installation and processor run procedures are different in case of developers and basic users. The differences are underlined when needed. 


Prerequisites
^^^^^^^^^^^^^

* Conda should be already installed
* In case you are a developer
    * `git <https://tortoisegit.org>`_ should be already installed
    * `tortoisegit <https://tortoisegit.org>`_ optional, a git GUI for Windows
    * A Python IDE (i.e. `spyder <www.spyder-ide.org>`_ ,  `vs code <https://code.visualstudio.com>`_)


BioPAL installation default option: "pip install" (users)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

User default installation method.
BioPAL will be automatically downoladed from `pypi`_.

.. _pypi: `https://pypi.org`

In a conda command window, type the following instruction, which creates an empty biopal environment with no packages but with the correct python version installed (you can customize the environment name, modifying the "biopal" string)::

    conda create --name biopal python==3.7.1
		
Before executing pip, install GDAL library with conda, by executing following commands in a conda command window (first activate the created environment, than install GDAL inside)::

    conda activate biopal
    conda install GDAL

Now the "biopal" environment is ready for installation; install the package by executing following command::	

    pip install biopal


BioPAL installation for developers only: "pip install -e"
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The code is editable, thanks to the "-e" option (so for users is suggested the default installation option).

Make a local clone forking the `repository <https://github.com/BioPAL/BioPAL>`_ from the web interface, executing the clone command in a conda command window (or use the tortoisegit GUI)::

    git clone --branch <branchname> <remote-repo-url>

where:

* <remote-repo-url> is https://github.com/your_name_here/BioPAL.git (write your specific name)
* <branchname> is the branch to be cloned: currently there is only a branch called `main` so the clone command will be (write your specific name)::

In a conda command window, type the following instruction, which creates an empty biopal environment with no packages but with the correct python version installed (you can customize the environment name, modifying the "biopal" string)::
	
    conda create --name biopal python==3.7.1
		
Before executing pip, install GDAL library with conda, by executing following commands in a conda command window (first activate the created environment, than install GDAL inside)::

    conda activate biopal
    conda install GDAL

Now the "biopal" environment is ready for installation; 
first enter inside the /BioPAL folder, 
than install the package by executing following command (the "." after "-e" option means "current folder")::

    pip install -e .


BioPAL datasets
^^^^^^^^^^^^^^^

BioPAL gives easy access to several datasets that are used for examples in the documentation and testing. These datasets are hosted on our FTP server and must be downloaded for use. Contact <biopal@esa.int> to receive access to the dataset and for more information.


Setup Configuration
^^^^^^^^^^^^^^^^^^^

Open the `inputs/Input_File.xml` and update following sections withj absolute paths:

* `dataset_query->L1C_repository`: this folder contains the data stacks to be processed
* `dataset_query->auxiliary_products_folder`: this folder contains auxiliary parameters related to the data stacks of the L1cRepository
* `output_specification->output_folder`: this is the folder where the output will be saved (each run corresponds to a sub-folder formatted with the current date time)

*NOTE: Sample data (L1C_repository) and auxiliaries (auxiliary_products_folder) can be obtained by writing to* <biopal@esa.int>.

GDAL paths configuration
""""""""""""""""""""""""

The BioPAL GDAL paths are automatically found by the processor after a correct installation procedure.
Also note that, under Windows this only works in CMD and not in PowerShell command window.
In case of problems or for particular user cases, it is possible to manually specify such paths, in this case edit `biopal/conf/Configuration_File.xml`, uncomment the "gdal" section and insert your absolute paths for:

* `gdal_path`: this is the folder containing the GDAL executables, usually in the `/bin` subfolder of GDAL environment (containing e.g., *gdalwarp*, *gdal_translate*,... )
* `gdal_enviroment_path`: this is the GDAL_DATA environment variable path

TIP: the above paths depend on your machine environment. GDAL has been automatically installed during the above procedure of conda environment creation; for a standard installation with conda, the paths should be found in paths similar to the following (where *xxx* is an alphanumeric string depending on the GDAL version installed)

Windows:

* gdal_path (i.e.): `C:\ProgramData\Anaconda3\pkgs\libgdal-xxx\Library\bin`
* gdal_enviroment_path (i.e.): `C:\ProgramData\Anaconda3\pkgs\libgdal-xxx\Library\share\gdal`

Linux:

* gdal_path (i.e.): `/home/user/.conda/envs/biopal/bin`
* gdal_enviroment_path (i.e.): `/home/user/.conda/pkgs/libgdal-xxx/share/gdal`


Run the processor
^^^^^^^^^^^^^^^^^

Set the `inputs/Input_File.xml` as desired, the `dataset_query` section is already filled with default L1C_date and geographic_boundaries_polygon, to be used with the DEMO DataSet from ESA.

Set the AGB, FH, FD, TOMO_FH configuration sections present in `biopal/conf/Configuration_File.xml` as desired (default configuration parameters alreasy present)

Than the procedure is different (developer, users), depending on the installation option used.

Run the processor for users
"""""""""""""""""""""""""""

In a conda command window, type the following instruction, which activates the biopal environment::
    
    conda activate biopal

In the same conda command window, from any folder, execute::

    biopal --conf conffolder inputfilexml

where:

* `inputfilexml`: path of the BioPAL xml input file (i.e. `/inputs` )
* `conffolder`:   path of the folder containing BioPAL xml configuration files (i.e. `biopal/conf/`)

With following command, default configurations are used::

    biopal inputfilexml

With following command, the biopal execution help will be shown::

    biopal

Run the processor for developers
""""""""""""""""""""""""""""""""

In a conda command window, type the following instruction, which activates the biopal environment::

    conda activate biopal

On the same conda command window execute::
        
    biopal --conf conffolder inputfilexml

where:
* `inputfilexml`: path of the BioPAL xml input file (i.e. `/inputs` )
* `conffolder`:   path of the folder containing BioPAL xml configuration files (i.e. `/biopal/conf/`)

With the following command, default configurations present in `biopal/conf/` are used::

    biopal inputfilexml
    
With the following command, the biopal execution help will be shown::

    biopal


Run the processor for developers with a script for debug
""""""""""""""""""""""""""""""""""""""""""""""""""""""""

To run the processor with a script for debug (developers only), and execute it in an IDE

Create a new *.py* file, with a text editor, with following content (where `yourPath/BioPAL` should be replaced with the folder where the BioPAL distribution has been git-cloned), and save it (i.e. `run_biopal_debug.py`)::

    from pathlib import Path
    import sys
    import os
    biopal_path = Path( 'yourPath/BioPAL' )
    sys.path.append( str(biopal_path) )
    os.chdir(biopal_path)
    from biopal.__main__ import biomassL2_processor_run
    input_file_xml_path = biopal_path.joinpath('inputs', 'Input_File.xml')
    conf_folder = biopal_path.joinpath( 'biopal','conf')
    biomassL2_processor_run(input_file_xml_path, conf_folder )

Read the Tutorla section for other scripts, and for manual execution of a BioPAL chain, step by step

Execute the `run_biopal_debug.py` script within your preferred IDE options (i.e. run, debug, breakpoints enabled....).
	
