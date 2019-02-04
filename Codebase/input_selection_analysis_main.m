%% Max Input Analysis
%%% By: Zach Gima 2018-11-19
%%% Main script to call input_selection_analysis for analyzing inputs selected for sensitivity analysis
clear all
close all
clc
%% Set Path to Pre-computed MAX Sensitivity Files and set output folder path
% Directory location for max sensitivity .mat files; *****make sure they only have
% the .mat files for the inputs in them

% Baseline A: Full Parameter Set (1 Group)
% senspath = '/Users/ztakeo/Documents/GitHub/OED_ParamID/Codebase/InputLibrary/MaxSensInputs/BaselineA_trim/Unformatted/';
% outputfolderpath  = 'Plots/InputAnalysis/BaselineA_trim/';

% Baseline B: Collinearity Only (1 Group)
senspath = '/Users/ztakeo/Documents/GitHub/OED_ParamID/Codebase/InputLibrary/MaxSensInputs/BaselineB_trim/Unformatted/';
outputfolderpath  = 'Plots/InputAnalysis/BaselineB_trim/';

% Baseline C: Collinearity + Sensitivity (2 Groups)    
% senspath = '/Users/ztakeo/Documents/GitHub/OED_ParamID/Codebase/InputLibrary/MaxSensInputs/OED_trim/Unformatted/';
% outputfolderpath  = 'Plots/InputAnalysis/Tmax60_trim/';

%%%%% Perturbation Case
%minus50
% senspath = '/Users/ztakeo/Documents/GitHub/OED_ParamID/Codebase/InputLibrary/MaxSensInputs/minus50/Unformatted/';
% outputfolderpath  = 'Plots/Perturb/InputAnalysis/minus50/';

%plus50
% senspath = '/Users/ztakeo/Documents/GitHub/OED_ParamID/Codebase/InputLibrary/MaxSensInputs/plus50/Unformatted/';
% outputfolderpath  = 'Plots/Perturb/InputAnalysis/plus50/';

if exist(outputfolderpath,'dir') == 7
    rmdir(outputfolderpath,'s'); % delete folder first (assuming it exists); this prevents the folder from keeping older .mat files
end
mkdir(outputfolderpath);

%% Load Directory of Sensitivity Files into Struct

sens_files = dir(senspath);
% on Mac, use line below to ignore '.','..', and '.DS_Store' files that
% are loaded into r 
sens_files=sens_files(~ismember({sens_files.name},{'.','..','.DS_Store'}));
num_exp = length(sens_files);

%% Call input_selection_analysis
input_selection_analysis(sens_files,senspath,outputfolderpath)