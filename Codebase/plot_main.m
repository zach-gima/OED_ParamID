%% Plot for perturbation analysis
close all
clear all
clc

% Load Resuts
result_folder = '/Users/ztakeo/Box Sync/HPC/HPC1/03-Dec-2018 20_08_19/';
result_filename = 'G2G1_nom.mat';
result_path = strcat(result_folder,result_filename);
load(result_path);

%Set new output folder (overwrite output folder saved in .mat file
% plotting on Mac, save plots to same folder as .mat result files
output_folder = result_folder;

% Load Parameters
run param/params_nominal
% For perturbation analysis specifically where perturbing away from nominal
truth_param = Nominal_param;

% % Select only the parameters that were identified
% theta_0_true = theta_0_true(sel_k);
% truth_param = truth_param(sel_k);
% ci95_full = ci95_full(sel_k);

% Call Plot function
Param_ID_plot(truth_param,theta_0_true,sel_k,paramID_out,ci95_full,t_paramID,rmse_final,output_folder,LM_options,bounds,alg_states)



