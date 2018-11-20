%% Max Input Analysis
%%% By: Zach Gima 2018-11-19
%%% Main script to call input_selection_analysis for analyzing inputs selected for sensitivity analysis
clear all
close all
clc
%% Set Path to Pre-computed MAX Sensitivity Files
% Directory location for max sensitivity .mat files; *****make sure they only have
% the .mat files for the inputs in them
% senspath = '/Users/ztakeo/Documents/GitHub/OED_ParamID/Codebase/InputLibrary/MaxSensInputs/Tmax45/Unformatted/';
senspath = '/Users/ztakeo/Documents/GitHub/OED_ParamID/Codebase/InputLibrary/MaxSensInputs/Tmax60/Unformatted/';


% Set Output folder path
% outputfolderpath  = 'Plots/InputAnalysis/Tmax45/';
outputfolderpath  = 'Plots/InputAnalysis/Tmax60/';

%% Load Directory of Sensitivity Files into Struct

sens_files = dir(senspath);
% on Mac, use line below to ignore '.','..', and '.DS_Store' files that
% are loaded into r 
sens_files=sens_files(~ismember({sens_files.name},{'.','..','.DS_Store'}));
num_exp = length(sens_files);

%% Call input_selection_analysis
input_selection_analysis(sens_files,senspath,outputfolderpath)