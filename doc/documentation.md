# BiomassL2 Prototype Processor, AGB only, 24 JULY 2020

The processor is based on a single xml main input file accessible to the user (conf/Input_File.xml), where only the AGB flow can be enabled. 
The AGB flow takes as inputs its specific xml configuration file, conf/ConfigurationFile_AGB .xml (editable by the user before launch) and a dedicated xml input file which instead is set up by the processor (based on the Main Input File conf/Input_File.xml).

Note: 
```sh
this demo contains only the AGB flow; any other flow, if enabled in the main input file, will be ignored, producing no output.
```

# AGB estimation

AGB estimation is based on following blocks:

**Stack based (the processing is done for each stack separately):**

  - Screen calibration / ground steering: applies phase calibration screens to the stack (where available) and does ground steering (terrain at zero elevation).
  - Ground cancellation: removes the contribution of the ground scattering from the stack.
  - Incidence angle calibration: removes the radiometric dependence on the local incidence angle.
  - Geocoding: projects on geographic map the processed data based on user inputs (Area & Dates, Product Sampling).

**Global (the processing combines the outputs from all the available stacks processed):**
  
  - ROI & CAL management + block management: defines the AGB processing grid, setting up blocks (areas where forest properties are assumed constant), ROIs (points where the model optimization is done) and initial CALs (points where the AGB is known a priori from external sources).
  - Block priority management: defines the processing order of overlapping blocks, starting from blocks with initial CALs.
  - Local AGB estimation: the estimation is done for each block, propagating results based on the order defined in the previous step.
  - Mosaicking: composes the complete map, merging all the estimated blocks.