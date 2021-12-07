%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    START     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Workspace
addpath(genpath(pwd))

% Figure Properties
set(0,'DefaultAxesFontName','Arial', 'Units', 'centimeters');
set(0,'DefaultTextFontName','Arial');
set(0,'DefaultFigurePaperUnits','centimeters');
set(0,'DefaultFigureUnits', 'centimeters');
set(0,'DefaultFigurePaperPositionMode','auto');
set(0,'DefaultFigurePosition', [5 5 16 8]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   ANALYSIS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% run the B vs VEG analysis
run read_modelOutput_B_VEG.m
gone

%% run the VEG vs VEG_nok2 analysis
run read_modelOutput_VEG_VEGnok2.m
gone

%% run the VEG vs VEG_noGW2Soil analysis
run read_modelOutput_VEG_VEGnoGW2Soil.m
gone

%% run the PFT vs B vs VEG analysis
run read_modelOutput_PFT.m
gone

%% run the parameter uncertainty (B, VEG & PFT)
run plot_paramUncertainty.m
gone

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      PREP      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% consistency check observations
run consistencyCheck_constraints.m
gone

%% plot study area & removed grids
run plot_studyArea
gone

%% plot PFT classes
run plot_PFTclasses.m
gone

%% plot Koeppen-Geiger zones
run plot_KoeppenGeiger.m
gone

%% plot parameter behavior
run plot_paramBehaviour.m
gone

