%% Plot for perturbation analysis
close all
clear all
clc

% Load Resuts
% result_folder = '/Users/ztakeo/Box Sync/HPC/HPC1/03-Dec-2018 20_08_19/';
result_folder = '/Users/ztakeo/Box Sync/HPC/HPC2/A-5/';
result_filename = 'full_nom.mat';
result_path = strcat(result_folder,result_filename);
load(result_path);

% %Set new output folder (overwrite output folder saved in .mat file
% % plotting on Mac, save plots to same folder as .mat result files
output_folder = result_folder;

% Load Parameters
run param/params_nominal
% For perturbation analysis specifically where perturbing away from nominal
truth_param = Nominal_param;

% % Select only the parameters that were identified
% theta_0_true = theta_0_true(sel_k); 
% truth_param = truth_param(sel_k);
% ci95_full = ci95_full(sel_k);

%% Supplement missing variables for simulations that didn't fully finish
run param/params_nominal
init_cond = 'nom';
theta_0 = Nominal_param;

%%% Baseline C Case:
% selection_vector = zeros(25,2);
% selection_vector(:,2) = [1;1;1;1;0;0;1;1;1;0;1;1;0;1;1;0;0;1;0;0;1;1;1;1;0];
% sel_k = find(selection_vector(:,2));

% Baseline A: All Case
selection_vector = zeros(25,1);
selection_vector(:,1) = [1;1;1;1;0;0;1;1;1;1;1;1;1;1;1;1;1;1;1;0;1;1;1;1;1];
sel_k = find(selection_vector(:,1));


% change this to indicate how many groups are being perturbed; e.g. for G2G1, num_groups = 2
num_perturbedgroups = 2;  % for running V_sim_debug or perturbation analysis

% Perturb parameters of interest all at the beginning
perturb_index = find(selection_vector(:,1)); % G1
% perturb_index = find(selection_vector(:,2)); %G1 & G2

perturb_factor = 0.95;
theta_0(perturb_index) = perturb_factor*theta_0(perturb_index);

theta_0_true = theta_0;


ci95_full = zeros(25,1);

t_paramID = [0;76734.75692511875];
rmse_final = [0.002309604133346; 1.000243447077960e-04];

% Set Levenberg-Marquardt Conditions
param_exit_thresh = 1e-6; % 1e-6 default values used in matlab (called StepTol)
chi_sq_exit_thresh = 1e-6; % 1e-6 default values used in matlab (called FuncTol)

LM_options.exit_cond = [param_exit_thresh, chi_sq_exit_thresh];
LM_options.maxIter = 20;
LM_options.ctrl_lambda =  [1e-2;1e-2];

alg_states = [];

run param/params_bounds % loads bounds struct (has fields bounds.min and bounds.max)


%% Call Plot function

Param_ID_plot(truth_param,theta_0_true,sel_k,paramID_out,ci95_full,t_paramID,rmse_final,output_folder,LM_options,bounds,alg_states,selection_vector)



