# VHFs [![DOI](https://zenodo.org/badge/638012149.svg)](https://zenodo.org/badge/latestdoi/638012149)
Dataset and Codes for Generating Figures in Aparco-Lara et al., 2023

![](https://img.shields.io/github/stars/JhonKevin3015/VHFs.svg) ![](https://img.shields.io/github/forks/JhonKevin3015/VHFs.svg) ![](https://img.shields.io/github/tag/JhonKevin3015/VHFs.svg) ![](https://img.shields.io/github/release/JhonKevin3015/VHFs.svg) ![](https://img.shields.io/github/issues/JhonKevin3015/VHFs.svg) 

## DATA

* The `VHFs/DATA` folder contains all the necessary data for the scripts found in the `MATLABcodes` and `PYTHONcodes` directories. The data is available in both MATLAB and Python formats.

* Inside the `DATA/Experiments/EXP*` directories, you will find the output data from experiments for both the unforced and forced cases on the 15th day of simulation.

* The `DATA/grid` directory contains the grid dimensions that were used in the simulation.

## MATLABcodes

* In `VHFs/DATA` folder are all data needed in the scripts contained in `MATLABcodes` and `PYTHONcodes`. The data is in both MATLAB and Python formats.
  - The `DATA/Experiments/EXP*` folders contain the outputs of experiments for the 15th day of simulation for both unforced and forced cases. And also contain the MITgcm numerical experiment configurations [more detail below](EXPE).
  - The `DATA/grid` folder contains the grid dimensions used in the simulation.

* This folder contains MATLAB codes to generate Figures 1, 2, 3, 4, 5, and 8.

* The `utils` folder includes [MITgcm tools](https://github.com/MITgcm/MITgcm/tree/master/utils) for reading the outputs of the [numerical experiments](https://github.com/JhonKevin3015/Experiments) found in the `DATA/EXP*` directories.

* All MATLAB code files have the following lines to ensure the use of utility packages:
  ```matlab
  addpath(genpath([dirPC '../MATLABcodes/utils/matlab']))
  addpath(genpath([dirPC '../MATLABcodes/utils/matlab_functions']))

(Note: Make sure to replace `dirPC` with the appropriate directory path before using the code.)


 ## PYTHONcodes

* This folder contains Python codes to generate Figures 6, 7, 9, 10, 11, and 12.

* The `SPECTRAL.py` and `wf_spectrum.py` files are Python functions.



# ------------------------------------------------------
<a id="EXPE"></a>
### DATA/Experiments

The MITgcm numerical experiment configurations utilized in the study by Aparco-Lara et al., 2023, are as follows:

####  EXP1 
Unforced case with KPP active
####	EXP2 
Forced case with KPP active
####	EXP3 
Unforced case with KPP inactive
####	EXP4 
Unforced case with KPP inactive
#### General Description
* Each experiment folder contains a `code` and `input` subfolder. Within the `EXP*/CODE` folder, you will find two important files:
  - `packages.conf`: This file determines whether KPP is activated or deactivated for the specific experiment. It is set to activate KPP for EXP1 and EXP2 and deactivate KPP for EXP3 and EXP4.
  - `SIZE.h`: This file sets the number of processors used for the simulations.


* In `EXP*/input` folder, you can find the data files to modify the physical parameters, as well as the location where the Initial Conditions files obtained from Matlab codes in the `IC` folder should be placed.

* Additionally, in the `EXP2` and `EXP4` folders, you will find an additional file named `QnetIDE.forcing`, which contains the atmospheric forcing data.

* Each experiment folder contains outputs for the 15th day of simulation. This data was used to [generate figures](https://github.com/JhonKevin3015/VHFs) in Aparco-Lara et al., 2023.

  
### IC: Matlab-code to Generate Initial Conditions and Qnet-Forcing Binary Files.

To set up the initial conditions and atmospheric forcing for the MITgcm numerical experiments, Matlab codes are utilized. Follow these steps to properly configure the necessary files:

* Run the Matlab script `InitCondi.m` to generate the binary file required for the initial conditions. This binary file should be placed in each `EXP*/input/` folder corresponding to the specific experiment (EXP1, EXP2, EXP3, EXP4).

* Additionally, the Matlab script `InitCondi.m` also generates the `QnetIDE.forcing` file. This file should be placed in the `EXP2/input/` and `EXP4/input/` folders, as they require an additional atmospheric forcing component.

By following these steps, you will have prepared the necessary input files for the MITgcm numerical experiments. For more detailed guidance on running the simulations, you can refer to the [MITgcm tutorial (examples)](https://mitgcm.readthedocs.io/en/latest/examples/examples.html). This tutorial provides comprehensive information on running MITgcm simulations and can be a valuable resource for your research.




