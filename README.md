# BioPAL

# test

[![Documentation Status](https://readthedocs.org/projects/biopal/badge/?version=latest)](http://biopal.readthedocs.io/?badge=latest)
[![PyPI](https://img.shields.io/pypi/v/biopal)](https://pypi.org/project/biopal)
[![PyPI - License](https://img.shields.io/pypi/l/biopal)](https://pypi.org/project/biopal)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/biopal)](https://pypi.org/project/biopal)


The BIOMASS Product Algorithm Laboratory hosts official tools for processing and analysing ESA\'s BIOMASS mission data.

-   Website: www.biopal.org
-   Documentation:
-   Mailing: <biopal@esa.int>
-   Contributing:
-   Bug reports:

# Objective

BIOMASS is ESA's (European Space Agency) seventh Earth Explorer mission, currently scheduled for launch in 2023. The satellite will be the first P-band SAR (Synthetic Aperture Radar) sensor in space and will be operated in fully polarimetric interferometric and tomographic modes. The mission main aim is to map forest properties globally, but the sensor will also allow exploring subsurface scenarios (ice, desert).

The BIOMASS Product Algorithm Laboratory (BioPAL) is an evolution of the software developed for the [BIOMASS prototype processor](https://www.mdpi.com/2072-4292/12/6/985) into an open source library to be used and contributed by the scientific community.

This repository collects the software routines for processing Level 1 SAR products to generate Level 2 forest products of Above Ground Biomass (AGB), Forest Heigth (FH) and Forest Disturbance (FD). More details about these products and BIOMASS can be found [here](https://www.mdpi.com/2072-4292/12/6/985).

# Structure of the Project

This repository is organized as follows:

-   **arepytools**: Aresys I/O library for reading and managing the input dataset. Will be turned to an independent library in the future.

-   **biopal**: contains the BioPAL source code in particular:

    -   the`biopal/_package_data` folder (do not edit) contains the default Input and Configuration xml files (use biopal-quickstart to get editable ones, see **Getting Started** section below)

-   **doc**: contains the documentation.

BioPAL is already used by some ESA sponsored project, however it is still an experimental code.
This means that there might be still bugs. If you happen to find one, make us happy by filling an issue with all the details to reproduce it the best you can.

You can follow the developement roadmap to version 1.0.0 [here](https://github.com/BioPAL/BioPAL/projects/2).


# Getting Started

For advanced insatallation and usage options, refer to **Documentation** section below.

## BioPAL installation (default option)
This installation procedure makes use of the open-source package management system [conda](https://docs.conda.io/projects/conda/en/latest/), to be pre-installed.

Open a command window with *conda* available and follow this procedure.

Create an empty biopal environment:

    conda create --name biopal python==3.7.1

Install GDAL library:

    conda activate biopal
    conda install GDAL

Install the package:

    pip install biopal

Configure biopal:

    biopal-quickstart FOLDER

"FOLDER" is the path where usable and editable versions of `Input_File.xml` and `Configuration_File.xml` files will be generated.

## Run BioPAL

Prepare your `Input_File.xml` and `Configuration_File.xml`, than open a command window with *conda* available and run BioPAL:

    conda activate biopal
    biopal --conf conf_folder inputfilexml

* *inputfilexml*: path of the `Input_File.xml` 
* *conf_folder*:  path of the folder containing `Configuration_File.xml`


# BioPAL datasets

BioPAL gives easy access to several datasets that are used for examples in the documentation and testing. 
These datasets are hosted on our FTP server and must be downloaded for use. 

Contact <biopal@esa.int> to receive access to the dataset and for more information.


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

Documentation is work in progress and can be found in  [project doc/ folder](https://github.com/BioPAL/BioPAL/tree/main/doc) and on [website doc section](https://www.biopal.org/docs/).

The user manual of the previous prototype software can be found in [legacy](https://github.com/BioPAL/BioPAL/tree/main/doc/legacy/ARE-017082_BIOMASS_L2_User_Manual_[prototype_legacy].pdf).

# History

BioPAL was originally written and is currently maintained by Aresys and the BioPAL team on behalf of ESA.

BioPAL team includes reperesentatives of several european research institutions, see [website about section](https://www.biopal.org/about/).


# Citing

If you use BioPAL, please add these citations:

-   *BioPAL: BIOMASS Product Algorithm Laboratory, https://github.com/BioPAL/BioPAL*

-   *Banda F, Giudici D, Le Toan T, Mariotti d’Alessandro M, Papathanassiou K, Quegan S, Riembauer G, Scipal K, Soja M, Tebaldini S, Ulander L, Villard L. The BIOMASS Level 2 Prototype Processor: Design and Experimental Results of Above-Ground Biomass Estimation. Remote Sensing. 2020; 12(6):985. https://doi.org/10.3390/rs12060985*