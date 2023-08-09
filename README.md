# VHFs [![DOI](https://zenodo.org/badge/638012149.svg)](https://zenodo.org/badge/latestdoi/638012149)
Dataset and Codes for Generating Figures in Aparco-Lara et al., 2023

![](https://img.shields.io/github/stars/JhonKevin3015/VHFs.svg) ![](https://img.shields.io/github/forks/JhonKevin3015/VHFs.svg) ![](https://img.shields.io/github/tag/JhonKevin3015/VHFs.svg) ![](https://img.shields.io/github/release/JhonKevin3015/VHFs.svg) ![](https://img.shields.io/github/issues/JhonKevin3015/VHFs.svg) 

### DATA

* The `VHFs/DATA` folder contains all the necessary data for the scripts found in the `MATLABcodes` and `PYTHONcodes` directories. The data is available in both MATLAB and Python formats.

* Inside the `DATA/EXP*` directories, you will find the output data from experiments for both the unforced and forced cases on the 15th day of simulation.

* The `DATA/grid` directory contains the grid dimensions that were used in the simulation.

### MATLABcodes

* In `VHFs/DATA` folder are all data needed in the scripts contained in `MATLABcodes` and `PYTHONcodes`. The data is in both MATLAB and Python formats.
  - The `DATA/EXP*` folders contain the outputs of experiments for the 15th day of simulation for both unforced and forced cases.
  - The `DATA/grid` folder contains the grid dimensions used in the simulation.

* This folder contains MATLAB codes to generate Figures 1, 2, 3, 4, 5, and 8.

* The `utils` folder includes [MITgcm tools](https://github.com/MITgcm/MITgcm/tree/master/utils) for reading the outputs of the [numerical experiments](https://github.com/JhonKevin3015/Experiments) found in the `DATA/EXP*` directories.

* All MATLAB code files have the following lines to ensure the use of utility packages:
  ```matlab
  addpath(genpath([dirPC '../MATLABcodes/utils/matlab']))
  addpath(genpath([dirPC '../MATLABcodes/utils/matlab_functions']))

(Note: Make sure to replace `dirPC` with the appropriate directory path before using the code.)


 ### PYTHONcodes

* This folder contains Python codes to generate Figures 6, 7, 9, 10, 11, and 12.

* The `SPECTRAL.py` and `wf_spectrum.py` files are Python functions.

