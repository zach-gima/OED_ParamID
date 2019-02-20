%% Simulate DFN Output Voltage
% Script to call & unit test DFN_sim_noRc_ZTG.m, sens_lm_cal..., and
% DFN_sim_casadi (ZTG 2018-6-11)

% Uncomment Group (1,2,3,4) you want to simulate voltage for
% If concatenating inputs, Just manually combine. e.g. G2G1 combine G2 & G1

% Also used for unit testing
%%%NOTE: DFN_sim_noRc must take in a SelecParam vector of 18x1 (i.e.
%%%excluding the equilibrium parameters). Otherwise, park0 will be assigned
%%%incorrectly (parameters will be assigned the valeus for other
%%%parameters, and IDAS will throw an error because  the integration will
%%%yield impossible values

% addpath('/Users/ztakeo/Documents/MATLAB/casadi')
% addpath('C:/Users/Zach/Documents/MATLAB/casadi_windows')
% import casadi.*

clearvars
close all
clc

%% Load Input
Inputs = load('InputLibrary/ValidationCycles/validation_cycles.mat');
% T_amb_sim = 25;
Current_exp = Inputs.Current_exp;
Time_exp = Inputs.Time_exp;
V_LM_CELL = Inputs.V_LM_CELL;
T_amb_sim = Inputs.T_amb_sim; % note comes in celcius
exp_num = Inputs.exp_num;

Num_inputs = length(V_LM_CELL);

% In experimental ID, Rc needs to be identified for each experiment
% (Rc_tune) and the Rc value must be attached to each profile
if isfield(Inputs,'Rc')
    Rc = Inputs.Rc;
else %M2M case
    Rc = cell(length(Current_exp),1); % Create a cell num_exp x 1
    Rc(:,1) = {p.R_c}; % For M2M case, just use nominal Rc value 
end
    
%% Load Parameters
run param/params_NCA % loads p struct

%%%%% For Generating Truth Data
% run param/params_truth % load true param values vector: truth_param
% theta_0 = truth_param;
% 
% % Replace nominal param values w/ true values (for generating truth data)
% p.D_s_n0 = truth_param(1);
% p.D_s_p0 = truth_param(2);
% p.R_s_n = truth_param(3);
% p.R_s_p = truth_param(4);
% p.epsilon_s_n = truth_param(5); %  Equil. Struct
% p.epsilon_s_p = truth_param(6); %  Equil. Struct
% p.sig_n = truth_param(7);
% p.sig_p = truth_param(8);
% p.ElecFactorD = truth_param(9);
% p.epsilon_e_n = truth_param(10);
% p.epsilon_e_s = truth_param(11);
% p.epsilon_e_p = truth_param(12);
% p.ElecFactorK = truth_param(13);
% p.t_plus = truth_param(14);
% p.ElecFactorDA = truth_param(15);
% p.k_n0 = truth_param(16);
% p.k_p0 = truth_param(17);
% p.R_f_n = truth_param(18);
% p.R_f_p = truth_param(19);
% p.n_Li_s = truth_param(20); %  Equil. Struct
% p.c_e0 = truth_param(21);


%%%%% For Generating Validation Cycle Voltage Response for model-to-model
%%%%% comparison
run param/params_expID % load true param values vector: truth_param
theta_0 = expID_param;
% 
% % Use nominal param values
% 
% % Overwrite the nominal values in the param struct w/ true values
% p.D_s_n0 = theta_0(1); %[G2]
% p.D_s_p0 = theta_0(2); %[G2]
% p.R_s_n = theta_0(3); %[G1]
% p.R_s_p = theta_0(4); %[G1]
% % p.epsilon_s_n = theta_0(5); %  Equil. Struct
% % p.epsilon_s_p = theta_0(6); %  Equil. Struct
% p.sig_n = theta_0(7);
% p.sig_p = theta_0(8);
% p.ElecFactorD = theta_0(9); %[G2]
% p.epsilon_e_n = theta_0(10); %[G2]
% p.epsilon_e_s = theta_0(11);
% p.epsilon_e_p = theta_0(12);
% p.ElecFactorK = theta_0(13); %[G2]
% p.t_plus = theta_0(14);
% p.ElecFactorDA = theta_0(15); %[G2]
% p.k_n0 = theta_0(16);
% p.k_p0 = theta_0(17);
% p.R_f_n = theta_0(18);
% p.R_f_p = theta_0(19);
% % p.n_Li_s = theta_0(20); %  Equil. Struct
% p.c_e0 = theta_0(21);
% p.E.Dsn = theta_0(22);
% p.E.Dsp = theta_0(23);
% p.E.kn = theta_0(24);
% p.E.kp = theta_0(25);
% 
% %%% Update Dependencies
% % Specific interfacial surface area
% p.a_s_n = 3*p.epsilon_s_n / p.R_s_n;  % Negative electrode [m^2/m^3]
% p.a_s_p = 3*p.epsilon_s_p / p.R_s_p;  % Positive electrode [m^2/m^3]
% 
% % make element to caclulate phi_{s} by Saehong Park 
% p.epsilon_f_n = 1 - p.epsilon_s_n - p.epsilon_e_n;  % Volume fraction of filler in neg. electrode
% p.epsilon_f_p = 1 - p.epsilon_s_p - p.epsilon_e_p;  % Volume fraction of filler in pos. electrode

%% Load Inputs
%%% NO NEED TO LOAD/INCLUDE RC VALUES FOR M2M COMPARISON; ONLY NECESSARY IN
%%% EXPERIMENT VS. MODEL SETTING

% 2018-7-26 Special Case -- Testing why Sens. Matrix singular at beginning
% Inputs = load('m2m_comp/simulated_inputs/V_sim_G2G1.mat');
% Inputs = load('m2m_comp/simulated_inputs/V_sim_G1.mat');
% 
% Current_exp = Inputs.Current_exp;
% Time_exp = Inputs.Time_exp;
% Voltage_exp = Inputs.V_LM_CELL;
% y_dat = cell2mat(Voltage_exp);
% Num_inputs = length(Voltage_exp);

%%%%Debug
% Num_inputs = 1;
% load('InputLibrary/OED_Inputs/OED_Test_339m.mat');
% Current_exp1 = Current_exp(1:10);
% Time_exp1 = Time_exp(1:10);
% Voltage_exp1 = Voltage_exp(1);
% Rc_exp1 = Rc;

%%%G1
% Num_inputs = 5;
% 
% load('InputLibrary/OED_Inputs/OED_Test_339m.mat');
% Current_exp1 = Current_exp;
% Time_exp1 = Time_exp;
% Voltage_exp1 = Voltage_exp(1);
% % Rc_exp1 = Rc;
% 
% load('InputLibrary/OED_Inputs/OED_Test_342.mat');
% Current_exp2 = Current_exp;
% Time_exp2 = Time_exp;
% Voltage_exp2 = Voltage_exp(1);
% % Rc_exp2 = Rc;
% 
% load('InputLibrary/OED_Inputs/OED_Test_345.mat');
% Current_exp3 = Current_exp;
% Time_exp3 = Time_exp;
% Voltage_exp3 = Voltage_exp(1);
% % Rc_exp3 = Rc;
% 
% load('InputLibrary/OED_Inputs/OED_Test_348.mat');
% Current_exp4 = Current_exp;
% Time_exp4 = Time_exp;
% Voltage_exp4 = Voltage_exp(1);
% % Rc_exp4 = Rc;
% 
% load('InputLibrary/OED_Inputs/OED_Test_365m.mat');
% Current_exp5 = Current_exp;
% Time_exp5 = Time_exp;
% Voltage_exp5 = Voltage_exp(1);
% % Rc_exp5 = Rc;

%%%G2 (339, 345, and 365 already loaded)
% Num_inputs = 2;
% 
% load('InputLibrary/OED_Inputs/OED_Test_90.mat');
% Current_exp1 = Current_exp;
% Time_exp1 = Time_exp;
% Voltage_exp1 = Voltage_exp(1);
% 
% load('InputLibrary/OED_Inputs/OED_Test_120.mat');
% Current_exp2 = Current_exp;
% Time_exp2 = Time_exp;
% Voltage_exp2 = Voltage_exp(1);

%%%G3 (365 already loaded)
% Num_inputs = 4;
% 
% load('InputLibrary/OED_Inputs/OED_Test_292m.mat');
% Current_exp1 = Current_exp;
% Time_exp1 = Time_exp;
% Voltage_exp1 = Voltage_exp(1);
% 
% load('InputLibrary/OED_Inputs/OED_Test_467.mat');
% Current_exp2 = Current_exp;
% Time_exp2 = Time_exp;
% Voltage_exp2 = Voltage_exp(1);
% 
% load('InputLibrary/OED_Inputs/OED_Test_503.mat');
% Current_exp3 = Current_exp;
% Time_exp3 = Time_exp;
% Voltage_exp3 = Voltage_exp(1);
% 
% load('InputLibrary/OED_Inputs/OED_Test_507.mat');
% Current_exp4 = Current_exp;
% Time_exp4 = Time_exp;
% Voltage_exp4 = Voltage_exp(1);
 
%%%G4(120, 503, 507 already loaded)
% Num_inputs = 2; 
% 
% load('InputLibrary/OED_Inputs/OED_Test_343.mat');
% Current_exp1 = Current_exp;
% Time_exp1 = Time_exp;
% Voltage_exp1 = Voltage_exp(1);
% 
% load('InputLibrary/OED_Inputs/OED_Test_362m.mat');
% Current_exp2 = Current_exp;
% Time_exp2 = Time_exp;
% Voltage_exp2 = Voltage_exp(1);

% clear Voltage_exp Temp_exp Description Current_exp Time_exp

%% Setup for Optimal experiments
% %DFN_sim_no_Rc works by taking in the sensitivity selection vector (SensSelec) and the
% %parameter vector (truth_param) to use. Based on which parameters are selected (here, all of
% %them), the function will replace the nominal value with the one specified
% %by the parameter vector passed into the function
% % selection_vector = [0;0;1;1;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0]; %G1
% % selection_vector = [1;1;1;1;0;0;0;0;1;1;0;0;1;0;1;0;0;0;0;0;0]; %G2
% % selection_vector = [1;1;1;1;0;0;0;0;1;1;0;1;1;0;1;1;0;1;1;0;1]; %G3
% selection_vector = [1;1;1;1;0;0;1;1;1;1;1;1;1;1;1;1;1;1;1;0;1]; %G4
% 
% % Selected_params = truth_param(sel_k);
% % S_LM_OLD = cell(Num_inputs,1);
% % S_LM_CASADI = cell(Num_inputs,1);
% 
% SensFlag = 0; %Turn sensitivity calculation off or on
% 
% % Create a composite cell of all the current, time, and voltage data for
% % each input used to identify the parameters in said group.
% 
% Current_exp = cell(Num_inputs,1);
% Time_exp = cell(Num_inputs,1);
% Voltage_exp = cell(Num_inputs,1);
% Rc_exp = cell(Num_inputs,1);
% 
% for i=1:Num_inputs 
%     Current_exp{i} = eval(['Current_exp' num2str(i)]);
%     Time_exp{i} = eval(['Time_exp' num2str(i)]);
%     Voltage_exp{i} = eval(['Voltage_exp' num2str(i)]);
% %     Rc_exp{i} = eval(['Rc_exp' num2str(i)]);
% end

%% Setup for Validation Cycles
selection_vector = [1;1;1;1;0;0;1;1;1;1;1;1;1;1;1;1;1;1;1;0;1;1;1;1;1]; %All, Year 2
SensFlag = 0; %Turn sensitivity calculation off or on
% 
% %%%%Validation Cycles
% Num_inputs = 16;
% exp_num = {'C1';'C2';'C3';'C4';'C5';'C6';'C7';'C8';'D1';'D2';'D3';'D4';'D5';'D6';'D7';'D8'};
% Current_exp_cell = cell(Num_inputs,1);
% Time_exp_cell = cell(Num_inputs,1);
% Voltage_exp_cell = cell(Num_inputs,1);
% T_amb_sim = cell(Num_inputs,1);
% 
% for ii = 1:Num_inputs
%    input_filename = strcat('InputLibrary/ValidationCycles/Unformatted/',exp_num{ii});
%    load(input_filename);
%    Current_exp_cell{ii} = Current_exp;
%    Time_exp_cell{ii} = Time_exp;
%    Voltage_exp_cell{ii} =  Voltage_exp; %for getting V0
%    T_amb_sim{ii} = 25;
% end
% 
% %Rename
% Current_exp = Current_exp_cell;
% Time_exp = Time_exp_cell;
% Voltage_exp = Voltage_exp_cell;
% 
% clear Current_exp_cell Time_exp_cell Voltage_exp_cell

%% Call DFN function to simulate voltage (& sensitivity)
% V_LM_CELL_sim = cell(Num_inputs,1);
% alg_states = cell(Num_inputs,1);
% cssn_sim = cell(Num_inputs,1);
% cssp_sim = cell(Num_inputs,1);
% etan_sim = cell(Num_inputs,1);
% etap_sim = cell(Num_inputs,1);
% T1_sim = cell(Num_inputs,1);
% T2_sim = cell(Num_inputs,1);

SensSelec = selection_vector;
sel_k = find(SensSelec);
Selected_params = theta_0(sel_k);

parfor idx = 1:Num_inputs
% for idx = 1:Num_inputs
%     V_LM_CELL{idx} =  DFN_sim_noRc_ZTG(Current_exp{idx}, Time_exp{idx}, Voltage_exp{idx},SensSelec,Selected_params); %, Rc_exp{idx});
%    V_true{idx} = DFN_sim_noRc(Current_exp{idx}, Time_exp{idx}, Voltage_exp{idx},SensSelec,Selected_params);

   [V_LM_CELL_sim{idx},~] = DFN_sim_casadi(p, exp_num{idx}, Current_exp{idx}, Time_exp{idx}, V_LM_CELL{idx}, T_amb_sim{idx}, SensSelec, Selected_params, SensFlag,Rc{idx}); % SensFlag == 0
   
   % Parse specific algebraic states of interest
%    cssn_sim{idx} = alg_states{idx}.cssn_sim;
%    cssp_sim{idx} = alg_states{idx}.cssp_sim;
%    etan_sim{idx} = alg_states{idx}.etan_sim;
%    etap_sim{idx} = alg_states{idx}.etap_sim;
%    T1_sim{idx} = alg_states{idx}.T1_sim;
%    T2_sim{idx} = alg_states{idx}.T2_sim;

   
%     [V_LM_CELL{idx}, S_LM_CASADI{idx}] = DFN_sim_casadi(p,Current_exp{idx}, Time_exp{idx}, Voltage_exp{idx}, SensSelec, Selected_params,SensFlag);
%     [V_sens{idx},S_LM_OLD{idx}] = sens_lm_cal_noRc_ZTG(Current_exp{idx}, Time_exp{idx}, Voltage_exp{idx},SensSelec,Selected_params);
end

% clear alg_states
% alg_states.cssn_sim = cssn_sim;
% alg_states.cssp_sim = cssp_sim;
% alg_states.etan_sim = etan_sim;
% alg_states.etap_sim = etap_sim;
save('validation_results.mat','V_LM_CELL_sim','Current_exp','Time_exp','V_LM_CELL','exp_num');

% save('V_sim_G2G1.mat','exp_num','Time_exp','T_amb_sim','T1_sim','T2_sim','Current_exp','V_LM_CELL','cssn_sim','cssp_sim','etan_sim','etap_sim');

% V_LM_CELL = cell2mat(V_LM_CELL);
% S_LM_CASADI = cell2mat(S_LM_CASADI);
% 
% V_sens = cell2mat(V_sens);
% S_LM_OLD = cell2mat(S_LM_OLD);

% fprintf('Determinant of sens_lm_cal_noRc_ZTG information matrix: %e \n',det(S_LM_OLD'*S_LM_OLD));
% fprintf('Determinant of DFN_sim_casadi information matrix: %e \n',det(S_LM_CASADI'*S_LM_CASADI));

% V_diff = V_sens - V_LM_CELL

%% uncomment to create .mat files where groups are combined G2G1, G3G2G1...
% S1 = load('V_sim_G1.mat');
% Current_exp_1 = S1.Current_exp;
% Time_exp_1 = S1.Time_exp;
% V_LM_CELL_1 = S1.V_LM_CELL;
% cssn_sim_1 = S1.cssn_sim;
% cssp_sim_1 = S1.cssp_sim;
% etan_sim_1 = S1.etan_sim;
% etap_sim_1 = S1.etap_sim;
% 
% S2 = load('V_sim_G2.mat');
% Current_exp_2 = S2.Current_exp;
% Time_exp_2 = S2.Time_exp;
% V_LM_CELL_2 = S2.V_LM_CELL;
% cssn_sim_2 = S2.cssn_sim;
% cssp_sim_2 = S2.cssp_sim;
% etan_sim_2 = S2.etan_sim;
% etap_sim_2 = S2.etap_sim;
% 
% S3 = load('V_sim_G3.mat');
% Current_exp_3 = S3.Current_exp;
% Time_exp_3 = S3.Time_exp;
% V_LM_CELL_3 = S3.V_LM_CELL;
% cssn_sim_3 = S3.cssn_sim;
% cssp_sim_3 = S3.cssp_sim;
% etan_sim_3 = S3.etan_sim;
% etap_sim_3 = S3.etap_sim;
% 
% S4 = load('V_sim_G4.mat');
% Current_exp_4 = S4.Current_exp;
% Time_exp_4 = S4.Time_exp;
% V_LM_CELL_4 = S4.V_LM_CELL;
% cssn_sim_4 = S4.cssn_sim;
% cssp_sim_4 = S4.cssp_sim;
% etan_sim_4 = S4.etan_sim;
% etap_sim_4 = S4.etap_sim;
% 
% filename = {'V_sim_G2G1.mat',;'V_sim_G3G2G1.mat';'V_sim_G4G3G2G1.mat'};
% 
% % G2G1
% Current_exp =  vertcat(Current_exp_2,Current_exp_1);
% Time_exp =  vertcat(Time_exp_2,Time_exp_1);
% V_LM_CELL = vertcat(V_LM_CELL_2,V_LM_CELL_1);
% cssn_sim = vertcat(cssn_sim_2,cssn_sim_1);
% cssp_sim = vertcat(cssp_sim_2,cssp_sim_1);
% etan_sim = vertcat(etan_sim_2,etan_sim_1);
% etap_sim = vertcat(etap_sim_2,etap_sim_1);
% 
% save(filename{1},'Current_exp','Time_exp','V_LM_CELL','cssn_sim','cssp_sim','etan_sim','etap_sim')
% clear Current_exp Time_exp V_LM_CELL 'cssn_sim' 'cssp_sim' 'etan_sim' 'etap_sim'
% 
% %G3G2G1
% Current_exp =  vertcat(Current_exp_3,Current_exp_2,Current_exp_1);
% Time_exp =  vertcat(Time_exp_3,Time_exp_2,Time_exp_1);
% V_LM_CELL = vertcat(V_LM_CELL_3,V_LM_CELL_2,V_LM_CELL_1);
% cssn_sim = vertcat(cssn_sim_3,cssn_sim_2,cssn_sim_1);
% cssp_sim = vertcat(cssp_sim_3,cssp_sim_2,cssp_sim_1);
% etan_sim = vertcat(etan_sim_3,etan_sim_2,etan_sim_1);
% etap_sim = vertcat(etap_sim_3,etap_sim_2,etap_sim_1);
% 
% save(filename{2},'Current_exp','Time_exp','V_LM_CELL','cssn_sim','cssp_sim','etan_sim','etap_sim')
% clear Current_exp Time_exp V_LM_CELL 'cssn_sim' 'cssp_sim' 'etan_sim' 'etap_sim'
% 
% %G4G3G2G1
% Current_exp =  vertcat(Current_exp_4,Current_exp_3,Current_exp_2,Current_exp_1);
% Time_exp =  vertcat(Time_exp_4,Time_exp_3,Time_exp_2,Time_exp_1);
% V_LM_CELL = vertcat(V_LM_CELL_4,V_LM_CELL_3,V_LM_CELL_2,V_LM_CELL_1);
% cssn_sim = vertcat(cssn_sim_4,cssn_sim_3,cssn_sim_2,cssn_sim_1);
% cssp_sim = vertcat(cssp_sim_4,cssp_sim_3,cssp_sim_2,cssp_sim_1);
% etan_sim = vertcat(etan_sim_4,etan_sim_3,etan_sim_2,etan_sim_1);
% etap_sim = vertcat(etap_sim_4,etap_sim_3,etap_sim_2,etap_sim_1);
% 
% save(filename{3},'Current_exp','Time_exp','V_LM_CELL','cssn_sim','cssp_sim','etan_sim','etap_sim')
% clear Current_exp Time_exp V_LM_CELL 'cssn_sim' 'cssp_sim' 'etan_sim' 'etap_sim'