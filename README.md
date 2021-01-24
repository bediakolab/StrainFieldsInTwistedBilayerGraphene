# StrainFieldsInTwistedBilayerGraphene
Matlab scripts and functions used in the publication "Strain Fields in Twisted Bilayer Graphene" (arXiv:2008.09761)

Bediako Lab, UC Berkeley, 01/23/2021

Contact information:
bediako@berkeley.edu (Prof. Kwabena Bediako),
nkazmier@caltech.edu (Nathanael Kazmierczak),
kazmierczak314@gmail.com (Nathanael Kazmierczak)

------------------------------------------------------------------

This repository contains the code that was used in all stages of the data analysis process. Only some of the code is still needed to process a 4D-STEM interferometry dataset on twisted bilayer graphene. The following functions, scripts, and folders represent the core functionality used in the final publication:

---The main data analysis file is "DSC_BlinkingFit_Code/FourDSTEM_Analysis_Engine.m". This file is > 11,000 lines long and contains the key functions for enabling most parts of the analysis. It defines a handle class for each dataset loaded in. Class methods are then invoked to update the state of the dataset analysis object, populating instance variables such as the displacement fields and strain maps. Throughout the repository, there are many scripts that set parameters for each individual datasets and invoke the appropriate methods. The analysis for each dataset was performed in stages, with the data analysis object saved in between the execution of each script. There are two main stages, detailed below.

---Stage 1 is exemplified by "MaddieStuff/revised2_script_for_data_fitting.m". This file converts the ~30GB raw .h5 dataset into the ~30MB reduced .mat data by defining masks for the graphene overlap regions, integrating the disks, performing background subtraction, running the displacement vector fits, etc. (See Supporting Information for details.) The 'saveObject()' function writes out an objectdata.mat data file containing the handle class instance with all instance variables. Stage 1 analysis typically takes 15-30 minutes or so, mostly due to the displacement vector fitting.

---Stage 2 performs the detailed strain and disorder analysis of each dataset. Common tasks for each dataset include (1) performing geometry fitting (Figure 2), (2) performing strain mapping (Figure 3), (3) performing pseudostacking analysis (Figure 4), and (4) performing heterostrain analysis (Figure 6). For (1), see "GeometryStatistics/AA_Gaussian_Circle_scripts/" for examples of the simple scripts used to make Figure 2e, and "GeometryStatistics/Session2Scripts/" for examples of the more involved scripts used to get the geometry registration information needed to make the the strain maps. For (2), see "StrainCalculationsAndAnnealing/Session2/" for phase unwrapping and "Filters/TGVstrainmappingscripts/" for strain mapping proper. These scripts make real space maps such as those seen in Figure 3a-d. For systematic studies of strain parameters as a function of twist angle (Figure 3f-h), see "StrainCalculationsAndAnnealing/StrainCalculationScripts/StrainCalcsPostColin/". For (3), see the folder "Pseudostacking/" (Figure 4a). For (4), see the folder "StrainCalculationsAndAnnealing/Uniaxial/" (Figure 6a).


Some notes:

-Absolute filepaths were used throughout. These will need to be reworked to point to the correct directory on a different machine. 

-Before using this code, raw .dm4 datasets were loaded into py4DSTEM (https://github.com/py4dstem/py4DSTEM) and saved as .h5 files for ease of loading into Matlab.

-Our workflow primarily used Matlab R2016b, and assumes the presence of the statistics, optimization, and parallel computing toolboxes. The code will run without the parallel computing toolbox, but the displacement fitting functions will be significantly slower.


