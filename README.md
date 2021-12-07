# VEGPaper
These scripts reproduce the analysis and figures of Trautmann, T., Koirala, S., Carvalhais, N., Güntner, A., &amp; Jung, M. (2021). The importance of vegetation to understand terrestrial water storage variations. Hydrology and Earth System Sciences Discussions, 1-31.

% All scripts are written and run with MATLAB R2019a

% To run the scripts:
- open MATLAB from the folder, e.g. by opening getStarted.mat
- getStarted.mat sets all paths and runs all other scripts
- 
- if paths are set correctly, all files in the main folder can be run independently

% download the data from zenodo and unzip them and add the folder data/ to the working directory
% data/ should contain:
  * data/input 	- forcing data and observational constraints 
  						  - ancillary data such as the regional classification
  						  - all data is provided in *.mat files* and adjusted to the study area and considered time period
  						  - all original data is publical available, please se the references and processing as described in Trautmann et al. 2021
  * data/model_output 	- output from different model experiments, including simulated fluxes & states and infomration on model settings
