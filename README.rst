BioPAL
======

The BIOMASS Product Algorithm Laboratory hosts official tools for processing and analysing ESA's BIOMASS mission data.

* Website:
* Documentation:
* Mailing: biopal@esa.int
* Contributing:
* Bug reports:

Objective
=========

BIOMASS is ESA’s (European Space Agency) seventh Earth Explorer mission, currently scheduled for launch in 2022.
The satellite will be the first P-band SAR (Synthetic Aperture Radar) sensor in space and will be operated in fully polarimetric interferometric and tomographic modes.
The mission main aim is to map forest properties globally, but the sensor will also allow exploring subsurface scenarios (ice, desert).

The BIOMASS Product Algorithm Laboratory (BioPAL) is an evolution of the software developed for the `BIOMASS prototype processor <https://www.mdpi.com/2072-4292/12/6/985>`_ into an open source library to be used and contributed by the scientific community.

This repository collects the software routines for processing Level 1 SAR products to generate Level 2 forest products of Above Ground Biomass (AGB). More details about this products and BIOMASS can be found `here <https://www.mdpi.com/2072-4292/12/6/985>`_.

Structure of the Project
========================

This repository is organized as follows:

* **arepytools**: I/O library for reading and managing the input dataset. Will be turned to an independent library in the future.
* **biopal**: contains the BioPAL source code in particular:

		* the `biopal/conf/ConfigurationFile_AGB.xml` contains the configuration for the Above Ground Biomass process.
		* the `biopal/conf/biopal_configuration.xml` contains the configuration for the BioPAL module environment.
		
* **doc**: contains the documentation.
* **inputs**: contains input data and the Input File, to be set by the user before running an instance of AGB processing.

Getting Started
===============

Here the required environment and dependencies are listed, as well as installation options available for now (installation options will be expanded including PyPi, Docker).

Requirements
------------

Python 3.7 is a minimum requirement, along with the following packages:

* Equi7Grid 0.0.10
* Pytileproj 0.0.12
* GDAL 2.3.3
* lxml 4.4.1
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
* scikit-image 0.16.2
* progressbar2 3.50.0

Installing without pip
----------------------

#. Download and unzip the current BioPAL distribution to your local hard drive.

#. Open a terminal and change directory to the folder that contains the setup.py file.

#. Install BioPAL by running::

    python setup.py install

Installing from Github
----------------------

Installation can be directly done from the master node (not recommended)::

	pip install git+https://git@github.com/BioPAL/BioPAL.git@master

Setup Configuration
-------------------

Open the `inputs/Input_File.xml` and update directory locations with absolute paths for:

* `L1cRepository`: this should correspond to the `inputs/dataSet` folder.
* `AuxiliaryProductsFolder`: this should correspond to the `inputs/auxiliary_data_pf` folder.
* `OutputFolder`: this is the folder where the output will be saved (each run corresponds to a sub-folder with the current date).

*NOTE: data and auxiliaries must be copied to the* \ `inputs` *folder. Sample data and auxiliaries can be obtained by writing to* biopal@esa.int.

Configure the BioPAL GDAL paths by editing `biopal/conf/biopal_configuration.xml` with absolute paths for:

* `gdal_path`: this is the folder containing the GDAL executables, usually in the `/bin` subfolder of GDAL environment (containing e.g., *gdalwarp*, *gdal_translate*, ... ).
* `gdal_enviroment_path`: this is the `GDAL_DATA` environment variable path.


BioPAL datasets
---------------

BioPAL gives easy access to several datasets that are used for examples in the documentation and testing. These datasets are hosted on our FTP server and must be downloaded for use. Contact biopal@esa.int to recieve access to the dataset and for more information.


Run
---

Be sure to configure the BioPAL GDAL paths as described above before.

#. Set as desired the `inputs/Input_File.xml`.
#. From terminal and with a Python environment installed, enter inside the main BioPAL folder and execute::
	
	python -m biopal
	
   this command automatically executes the `biopal/__main__.py` module, parsing the `inputs/Input_File.xml`.

Call for Contributions
======================

BioPAL is an open source project supported by a community who apprechiates help from a wide range of different backgrounds. Large or small, any contribution makes a big difference; and if you've never contributed to an open source project before, we hope you will start with BioPAL!

If you are interested in contributing, check out our contributor's guide. Beyond enhancing the processing algorithms, there are many ways to contribute:

   * Submit a bug report or feature request on GitHub Issues.
   * Contribute a Jupyter notebook to our examples gallery.
   * Assist us with user testing.
   * Add to the documentation or help with our website.
   * Write unit or integration tests for our project.    
   * Answer questions on our issues, slack channel, MAAP Forums, and elsewhere.
   * Write a blog post, tweet, or share our project with others.
   * Teach someone how to use BioPAL.

As you can see, there are lots of ways to get involved and we would be very happy for you to join us! The only thing we ask is that you abide by the principles of openness, respect, and consideration of others as described in our Code of Conduct.

Contributing Guidelines in Brief
--------------------------------

Read carefully also contributor's guides before getting started.

#. Fork the repository.

#. Clone the private fork locally (execute the following command in your terminal)::

	git clone https://github.com/your_name_here/BioPAL/BioPAL.git

#. Follow the instructions specified in the documentation, make a demo run and compare with reference output. Make sure all tests are passed.

#. Add the main repository to the list of your remotes (in order to pull the latest changes before making local changes)::

	git remote add upstream https://github.com/BioPAL/BioPAL

#. Create a branch for local development.

#. Commit local changes and push local branch to the GitHub private fork.

#. Submit a pull request through the GitHub website to the `dev` branch of the main repository.

Pull Request Requirements
-------------------------

#. Include new tests for all the new routines developed.
#. Documentation should be updated accordingly.
#. Updated code must pass all the tests.

Documentation
=============

For more details on processing configuration, refer to the `user manual <https://github.com/BioPAL/BioPAL/doc/ARE-017082_BIOMASS_L2_User_Manual.pdf>`_  version 1.1 .

History
=======

BioPAL was originally written and is currently maintained by Aresys and the BioPAL team on behalf of ESA.

Citing
======

If you use BioPAL, please add a citation:

* *BioPAL: BIOMASS Product Algorithm Laboratory, https://github.com/BioPAL/BioPAL*

Contribution
============

* Francesco Banda
* Emanuele Giorgi
* Maciej Soja
* Stefanie Lumnitz
* Paolo Mazzuchelli
* Klaus Scipal
* Davide Giudici
* Clément Albinet

Affilliations
=============

TBD