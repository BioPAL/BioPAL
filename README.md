# BioPAL

The BIOMASS Product Algorithm Laboratory hosts official tools for processing and analysing ESA\'s BIOMASS mission data.

-   Website: www.biopal.org
-   Documentation:
-   Mailing: <biopal@esa.int>
-   Contributing:
-   Bug reports:

# Objective

BIOMASS is ESA's (European Space Agency) seventh Earth Explorer mission, currently scheduled for launch in 2023. The satellite will be the first P-band SAR (Synthetic Aperture Radar) sensor in space and will be operated in fully polarimetric interferometric and tomographic modes. The mission main aim is to map forest properties globally, but the sensor will also allow exploring subsurface scenarios (ice, desert).

The BIOMASS Product Algorithm Laboratory (BioPAL) is an evolution of the software developed for the [BIOMASS prototype processor](https://www.mdpi.com/2072-4292/12/6/985) into an open source library to be used and contributed by the scientific community.

This repository collects the software routines for processing Level 1 SAR products to generate Level 2 forest products of Above Ground Biomass (AGB) and Forest Heigth (FH). More details about these products and BIOMASS can be found [here](https://www.mdpi.com/2072-4292/12/6/985).

# Structure of the Project

This repository is organized as follows:

-   **arepytools**: Aresys I/O library for reading and managing the input dataset. Will be turned to an independent library in the future.

-   **biopal**: contains the BioPAL source code in particular:

    -   the`biopal/_package_data` folder (ndo not edit) contains the default Input and Configuration xml files (use biopal-quickstart to get editable ones, see after)
-   **doc**: contains the documentation.

BioPAL is already used by some ESA sponsored project, however it is still an experimental code.
This means that there might be still bugs. If you happen to find one, make us happy by filling an issue with all the details to reproduce it the best you can.

You can follow the developement roadmap to version 1.0.0 [here](https://github.com/BioPAL/BioPAL/projects/2).

# Getting Started

Here the required environment and dependencies are listed, as well as installation options available for now.

## Requirements

Python >= 3.7.1
The needed python packages are specified in the file [requirements.txt](https://github.com/BioPAL/BioPAL/blob/main/requirements.txt).

## Installation

Installation procedure described here makes use of the open-source package management system [conda](https://docs.conda.io/projects/conda/en/latest/). 

Note that the installation and processor run procedures are different in case of developers and basic users. The differences are underlined when needed. 

##### Prerequisites

- Conda should be already installed
- In case you are a developer
  - [git](https://git-scm.com/downloads) should be already installed
  - [tortoisegit](https://tortoisegit.org/), optional, a git GUI for Windows
  - Python IDE (i.e. [spyder](https://www.spyder-ide.org/) or [vs code](https://code.visualstudio.com/))


##### BioPAL installation default option: "pip install" (users)

BioPAL will be automatically downoladed from [pypi](`https://pypi.org/`)

Open a command window with *conda* available and follow this procedure.

Create an empty biopal environment: you can customize the *biopal* environment name:
	
        conda create --name biopal python==3.7.1
		
Install GDAL library:

        conda activate biopal
        conda install GDAL

Install the package:		

        pip install biopal


##### BioPAL installation for developers only: "pip install -e"

The code is editable, thanks to the "-e" option.


First, you need to fork your own copy of BioPAL from [github BioPAL](https://github.com/BioPAL/BioPAL). Your private fork url will be something like https://github.com/your_name_here/BioPAL

Than open a command window with conda available and follow this procedure.

Make a local clone inside an empty folder (your_installation_folder/):

        cd your_installation_folder/
        git clone --branch <branchname> <remote-repo-url>

   - `<remote-repo-url>` is the url created during the fork from [github BioPAL](https://github.com/BioPAL/BioPAL).
   - `<branchname>` is the branch to be cloned: if it's a new fork, there is only a branch called *main*


Create an empty biopal environment: you can customize the biopal environment name:

        conda create --name biopal python==3.7.1

Install GDAL library:

        conda activate biopal
        conda install GDAL

Install the package:

        cd your_installation_folder/BioPAL/
        pip install -e .


## Setup Configuration

### Quick Start
Execute following command once to configure biopal:

        biopal-quickstart FOLDER

Where FOLDER is the path where editable `Input_File.xml` and `Configuration_File.xml` files will be placed.

### Edit Input and Configuration
Open the `Input_File.xml` and update following sections with absolute paths:

-   `dataset_query->L1C_repository`: this folder contains the data stacks to be processed
-   `dataset_query->auxiliary_products_folder`: this folder contains auxiliary parameters related to the data stacks of the L1cRepository
-   `output_specification->output_folder`: this is the folder where the output will be saved (each run corresponds to a sub-folder formatted with the current date time)

*NOTE: Sample data (L1C_repository) and auxiliaries (auxiliary_products_folder) can be obtained by writing to* <biopal@esa.int>.

#### GDAL paths configuration
The BioPAL GDAL paths are automatically found by the processor after a correct installation procedure.
Also note that, under Windows this only works in CMD and not in PowerShell command window.
In case of problems or for particular user cases, it is possible to manually specify such paths, in this case edit `biopal/Configuration_File.xml`, uncomment the "gdal" section and insert your absolute paths for:
-   `gdal_path`: this is the folder containing the GDAL executables, usually in the `/bin` subfolder of GDAL environment (containing e.g., *gdalwarp*, *gdal_translate*,... )
-   `gdal_enviroment_path`: this is the GDAL_DATA environment variable path

TIP: the above paths depend on your machine environment. GDAL has been automatically installed during the above procedure of conda environment creation; for a standard installation with conda, the paths should be found in paths similar to the following (where *xxx* is an alphanumeric string depending on the GDAL version installed)

###### Windows
- gdal_path (i.e.): `C:\ProgramData\Anaconda3\pkgs\libgdal-xxx\Library\bin`
- gdal_enviroment_path (i.e.): `C:\ProgramData\Anaconda3\pkgs\libgdal-xxx\Library\share\gdal`
###### Linux
- gdal_path (i.e.): `/home/user/.conda/envs/biopal/bin`
- gdal_enviroment_path (i.e.): `/home/user/.conda/pkgs/libgdal-xxx/share/gdal`
   

## BioPAL datasets

BioPAL gives easy access to several datasets that are used for examples in the documentation and testing. These datasets are hosted on our FTP server and must be downloaded for use. Contact <biopal@esa.int> to receive access to the dataset and for more information.

## Run the processor
Set the `Input_File.xml` as desired, the `dataset_query` section is already filled with default L1C_date and geographic_boundaries_polygon, to be used with the DEMO DataSet from ESA.

Set the AGB, FH, FD, TOMO_FH configuration sections present in `Configuration_File.xml` as desired (default configuration parameters alreasy present)

### Run the processor for users

Open a command window with conda available and follow this procedure.

Activate the biopal environment:

        conda activate biopal

Run BioPAL:

        biopal --conf conffolder inputfilexml

- `inputfilexml`: path of the BioPAL xml input file
- `conffolder`:   path of the folder containing BioPAL xml configuration file

Or Run BioPAL with default configurations:

        biopal inputfilexml

Or show BioPAL help:

        biopal -h

### Run the processor for developers with a script for debug

How to run the processor with a script to be launced from an IDE.

Create a new .py script file as (update needed paths):

        from pathlib import Path
        import sys
        import os
        biopal_path = Path( 'yourPath/BioPAL' )
        sys.path.append( str(biopal_path) )
        os.chdir(biopal_path)
        from biopal.__main__ import biomassL2_processor_run
        input_file_xml_path = biopal_path.joinpath('yourInputsPath', 'Input_File.xml')
        biomassL2_processor_run(input_file_xml_path, conf_folder='yourConfPath' )

Execute the script within your preferred IDE options (i.e. run, debug, breakpoints enabled....).

	
Read BioPAL [tutorial](https://www.biopal.org/docs/tutorials/biopal_first_tutorial/) for other examples to insert in the script


##### How to generate a Wheel package for pypi
1.  Verify not to have a ".nox" folder in BioPAL/.nox, otherwise delete it.

2.  From command window, with biopal environment enabled, digit:

    pip install nox                          (if not already installed in this environment)
    nox -s build_wheel -fb conda             (generate the wheel package; -fb conda needed only if anacondaconda is used)
    
	nox -s build_sdist -fb conda             (generate the sdist package, if needed)
	
	The build_wheel command will produce the wheel file "BioPAL\dist\biopal-0.1-py3-none-any.whl"
	The build_sdist command will produce the sdist (Software distribution)  file "BioPAL\dist\biopal-0.1.tar.gz"
    
	For pypi packages info see also https://python-packaging-tutorial.readthedocs.io/en/latest/uploading_pypi.html


# Call for Contributions

BioPAL is an open source project supported by a community who appreciates help from a wide range of different backgrounds. Large or small, any contribution makes a big difference; and if you\'ve never contributed to an open source project before, we hope you will start with BioPAL!

If you are interested in contributing, check out our contributor\'s guide. Beyond enhancing the processing algorithms, there are many ways to contribute:

-   Submit a bug report or feature request on GitHub Issues.
-   Contribute a Jupyter notebook to our examples gallery.
-   Assist us with user testing.
-   Add to the documentation or help with our website.
-   Write unit or integration tests for our project.
-   Answer questions on our issues, slack channel, MAAP Forums, and elsewhere.
-   Write a blog post, tweet, or share our project with others.
-   Teach someone how to use BioPAL.

As you can see, there are lots of ways to get involved and we would be very happy for you to join us! The only thing we ask is that you abide by the principles of openness, respect, and consideration of others as described in our Code of Conduct.

## Contributing Guidelines in Brief

Read carefully also contributor\'s guides before getting started.

1.  Fork the repository.

2.  Clone the private fork locally (execute the following command in your terminal):

        git clone https://github.com/your_name_here/BioPAL.git

3.  Follow the instructions specified in the documentation, make a demo run and compare with reference output. Make sure all tests are passed.

4.  Add the main repository to the list of your remotes (in order to pull the latest changes before making local changes):

        git remote add upstream https://github.com/BioPAL/BioPAL

5.  Create a branch for local development.

6.  Commit local changes and push local branch to the GitHub private fork.

7.  Submit a pull request through the GitHub website to the [main branch](https://github.com/BioPAL/BioPAL/tree/main) of the main repository.

## Pull Request Requirements

1.  Include new tests for all the new routines developed.
2.  Documentation should be updated accordingly.
3.  Updated code must pass all the tests.

# Documentation

Documentation is work in progress and can be found [here](https://github.com/BioPAL/BioPAL/tree/main/doc).

The user manual of the previous prototype software can be found in [legacy](https://github.com/BioPAL/BioPAL/tree/main/doc/legacy/ARE-017082_BIOMASS_L2_User_Manual_[prototype_legacy].pdf).

# History

BioPAL was originally written and is currently maintained by Aresys and the BioPAL team on behalf of ESA.

# Citing

If you use BioPAL, please add these citations:

-   *BioPAL: BIOMASS Product Algorithm Laboratory, https://github.com/BioPAL/BioPAL*

-   *Banda F, Giudici D, Le Toan T, Mariotti d’Alessandro M, Papathanassiou K, Quegan S, Riembauer G, Scipal K, Soja M, Tebaldini S, Ulander L, Villard L. The BIOMASS Level 2 Prototype Processor: Design and Experimental Results of Above-Ground Biomass Estimation. Remote Sensing. 2020; 12(6):985. https://doi.org/10.3390/rs12060985*

# Affilliations

TBD
