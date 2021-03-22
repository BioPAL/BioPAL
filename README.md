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

    -   the `biopal/conf/biopal_configuration.xml` contains the configuration for the BioPAL module environment;
    -   the`biopal/conf/ConfigurationFile_AGB.xml`contains the configuration for the Above Ground Biomass algorithm;
    -   the `biopal/conf/ConfigurationFile_FH.xml` contains the configuration for the interferometric phase Forest Heigth algorithm;
    -   the `biopal/conf/ConfigurationFile_TOMO_FH.xml` contains the configuration for the tomographic phase Forest Heigth algorithm.

-   **doc**: contains the documentation.

-   **inputs**: contains the XML Input File, to be set by the user before running an instance of the processing.

# Getting Started

Here the required environment and dependencies are listed, as well as installation options available for now.

## Requirements

Python 3.6 is a minimum requirement. The packages required are specified in the file [requirements.txt](https://github.com/BioPAL/BioPAL/blob/main/requirements.txt).

## Installation

Installation procedure described here makes use of the open-source package management system [conda](https://docs.conda.io/projects/conda/en/latest/). Note that the installation and processor run procedures are different in case of developers and basic users. The differences are underlined when needed. 

##### Prerequisites

- Conda should be already installed
- Python 3.6 or above should be already installed
- In case you are a developer
  - [git](https://git-scm.com/downloads) should be already installed
  - [tortoisegit](https://tortoisegit.org/), optional, a git GUI for Windows
  - Python IDE (i.e. [spyder](https://www.spyder-ide.org/)) to be installed in the environment (see next section for details)

##### Make a local clone (developers)

1. Fork the [repository](https://github.com/BioPAL/BioPAL) from the web interface.
2. Clone the private fork locally by executing the clone command in a conda command window (or use the tortoisegit GUI):

        git clone --branch <branchname> <remote-repo-url>
    where:
-- `<remote-repo-url>` is https://github.com/your_name_here/BioPAL.git (write your specific name)
-- `<branchname>` is the branch to be cloned: currently there is only a branch called `main`
    so the clone command will be (write your specific name):
    
        git clone --branch main https://github.com/your_name_here/BioPAL.git

##### Make a local clone (users)

1. If you are a user, download and unzip the current [BioPAL distribution](https://github.com/BioPAL/BioPAL) to your local hard drive.

##### Python environment preparation and activation

1.  In a conda command window, type the following instruction, which creates a biopal environment ( `environment.yml` is present into the BioPAL distribution ):
	
        conda env create --file environment.yml

2.  In the same conda command window, type the following instruction, which activates the created biopal environment:
	
        conda activate biopal

3.  *Only for developers*: In the same conda command window, install your favorite python IDE with a command as (in the following, *spyder* IDE is an example):
	
        conda install spyder

##### BioPAL package installation into the created environment (users)

If you are a *developer* you should not install the package: skip this part. You should refer to the section [Run the processor for developers](#-run-the-processor-for-developers). If you are a basic *user* follow this procedure. You should refer to the section [Run the processor for users](#-run-the-processor-for-users).

*NOTE: the installation procedure here described is not yet working, work in progress.*

In a conda command window, be sure to have the biopal environment active, otherwise digit

        conda activate biopal
        
then type the following instruction:

        pip install biopal

## Setup Configuration
Open the `inputs/Input_File.xml` and update directory locations with absolute paths for:

-   `L1cRepository`: this folder contains the data stacks to be processed
-   `AuxiliaryProductsFolder`: this folder contains auxiliary parameters related to the data stacks of the L1cRepository
-   `OutputFolder`: this is the folder where the output will be saved (each run corresponds to a sub-folder formatted with the current date time)

*NOTE: Sample data and auxiliaries can be obtained by writing to* <biopal@esa.int>.

The BioPAL GDAL paths are automatically found by the processor after a correct installation procedure; in case of problems or for particular user cases, it is possible to specify manually such paths, in this case edit `biopal/conf/biopal_configuration.xml`, uncomment the paths section and insert your absolute paths for:
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

## Run the processor for users
*NOTE: The run procedure here described is ready but not usable yet: it will be available once installation procedure described in "BioPAL package installation into the created environment" will be ready*.

1.  Set the `inputs/Input_File.xml` as desired
2.  Set the AGB, FH, TOMO_FH xml configurations present in `biopal/conf folder` as desired
3.  In a conda command window, type the following instruction, which activates the biopal environment:

        conda activate biopal

4.  In the same conda command window execute:

        biopal --conf conffolder inputfilexml
    where:
-- `inputfilexml`: path of the BioPAL xml input file (i.e. `/inputs` )
-- `conffolder`:   path of the folder containing BioPAL xml configuration files (i.e. `biopal/conf/`)
With following command, default configurations are used:

        biopal inputfilexml
    With following command, the biopal execution help will be shown:

        biopal

## Run the processor for developers
1.  Set the `inputs/Input_File.xml`as desired
2.  Set the AGB, FH, TOMO_FH xml configurations present in `biopal/conf` folder as desired
3.  In a conda command window, type the following instruction, which activates the biopal environment:

        conda activate biopal

    Then there are the following two choices: comand window or IDE.
##### To run the processor from command window (developers only):

4.  On the same conda command window enter in the BioPAL folder and execute:

        python -m biopal --conf conffolder inputfilexml
    where:
    -- `inputfilexml`: path of the BioPAL xml input file (i.e. `/inputs` )
    -- `conffolder`:   path of the folder containing BioPAL xml configuration files (i.e. `/biopal/conf/`)
    With the following command, default configurations present in `biopal/conf/` are used:

        python -m biopal inputfilexml
    With the following command, the biopal execution help will be shown:

        python -m biopal

##### To run the processor from a python IDE (developers only):
For this run method, a python IDE should be installed in the environment (see section [Python environment preparation and activation](#-python-environment-preparation-and-activation) above).

4. Create a new *.py* file, with a text editor, with following content (where `yourPath/BioPAL` should be replaced with the folder where the BioPAL distribution has been git-cloned), and save it (i.e. `run_biopal_debug.py`):

        from pathlib import Path
        import sys
        import os
        biopal_path = Path( `yourPath/BioPAL' )
        sys.path.append( str(biopal_path) )
        os.chdir(biopal_path)
        from biopal.__main__ import biomassL2_processor_run
        input_file_xml_path = biopal_path.joinpath('inputs', 'Input_File.xml')
        conf_folder = biopal_path.joinpath( 'biopal','conf')
        biomassL2_processor_run(input_file_xml_path, conf_folder )

5.  In the same conda command window open the python IDE (i.e. *spyder*) and execute the `run_biopal_debug.py` script with your preferred IDE options (i.e. run, debug, breakpoints enabled....).

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

For more details on processing configuration, refer to the initial release [user manual](https://github.com/BioPAL/BioPAL/doc/ARE-017082_BIOMASS_L2_User_Manual.pdf) version 1.1.

# History

BioPAL was originally written and is currently maintained by Aresys and the BioPAL team on behalf of ESA.

# Citing

If you use BioPAL, please add a citation:

-   *BioPAL: BIOMASS Product Algorithm Laboratory, https://github.com/BioPAL/BioPAL*

# Affilliations

TBD