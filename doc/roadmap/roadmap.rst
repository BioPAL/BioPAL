==============
BioPAL Roadmap
==============

This roadmap outlines tasks and features we will be investing resources
in to reach BioPAL's first releas, v0.1.0, in January 2021. It may be used
to encourage and inspire developers and collaboration.

Discussions about what will be included in the roadmap can be found:

* in `meeting_notes`_ 
* in the BioPAL  `roadmap_issue`_


Vision 2021
===========

- BioPAL strives to provide an example of Open Source Project best practices for the BIOMASS community, including governance structures, contributing guidlines, code contribution agreements, continuous integration and testing, and other collaborative infrastructure.
- BioPAL will contain a first open source version (v0.1.0) of the working BioL2 AGB processor, including improved documentation, code formatting and a user friendly Python API.
- The BioPAL community aims to refactor the current API into a higher and lower level user API, leveraging python modules, classes and functions. This API supports scientists to focus and improve on one aspect of the source code, through lower level access while full orchestration capabilities remain in place.


Tasks
=====

Community guidlines and Documentation
-------------------------------------


Checklist
~~~~~~~~~
- [ ] Governance structure
    - [x] document draft
    - [ ] legal check
- [ ] Code of conduct
    - [x] document
    - [x] committee
    - [ ] contact form
- [ ] Contribution guidlines
    - [ ] document for BioPAL
    - [ ] code contribution agreement
- [ ] Logo
    - [x] draft
    - [ ] final
    - [ ] stickers
- [ ] Website
    - [ ] domain
    - [x] first content draft
    - [ ] infrastructure
- [x] E-mail
- [x] Readme


API refactoring
---------------

* First step: create two main APPs, a first APP for stack based pre-processing and a second APP for global AGB estimation. This should be done by providing I/O to disk and configuration for each APP, as well as separating configuration from the functional execution in the code (classes design).
* Second step: create atomic APPs inside the two main ones. Suitable functionalities for encapsulation to be identified from the code.

Optimization
------------

* Improve Windows experience by shortening all the output and temporary folders/subfolders path length (max length for Windows is 250 characters).
* Update `get_min_time_stamp` function to get the time stamps from binary headers instead of xml metadata, this will speed-up the process.

Maintenance
-----------

* Analyze the warnings returned by the processor, improve the logging.


.. _`roadmap issue`: https://github.com/BioPAL/BioPAL/issues/2
.. _`meeting notes`: https://github.com/BioPAL/community/tree/master/00_dev_meetings