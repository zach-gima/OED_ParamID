clear all 
close all
clc

%% Load inputs
% Load Nominal Parameter Set, Bounds, & True params
run param/params_NCA % loads p struct
run param/params_bounds % loads bounds struct (has fields bounds.min and bounds.max)
run param/params_truth % loads truth_param array

filename_input = '/Users/ztakeo/Documents/GitHub/OED_ParamID/Codebase/InputLibrary/Experimental/V_sim_G1.mat';
Inputs = load(filename_input);

Current_exp = Inputs.Current_exp;
Time_exp = Inputs.Time_exp;
Voltage_exp = Inputs.V_LM_CELL;
T_amb = Inputs.T_amb_sim; % note comes in celcius
exp_num = Inputs.exp_num;
Rc = Inputs.Rc;

num_exp = length(Current_exp);
%%   Parameter Initial Conditions & Parameters to Identify
run param/params_nominal
init_cond = 'nom';
theta_0 = Nominal_param;

% Sets the index of the parameters to identify
selection_vector = [1;1;1;1;0;0;0;0;1;0;0;0;0;0;1;0;0;0;0;0;0;0;0;0;0]; %G1
SensSelec = selection_vector;
sel_k = find(selection_vector);
SelecParam = theta_0(sel_k);

%% Call DFN_sim_casadi
SensFlag = 0; % 0 = no sensitivity calculation, just voltage simulation

parfor idx = 1:length(num_exp)
    [V,~,~] = DFN_sim_casadi(p, exp_num{idx}, Current_exp{idx}, Time_exp{idx}, Voltage_exp{idx}, T_amb{idx}, SensSelec, SelecParam, SensFlag,Rc{idx}); % [ZTG change] removed Rc for no model-to-model comparison
end