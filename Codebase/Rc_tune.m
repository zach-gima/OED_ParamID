%%% Script for tuning the IR drop by setting Rc. This is necessary for
%%% identifying parameters from experimental data
clearvars
close all
clc
%% Load Parameters

run param/params_nominal
theta_0 = Nominal_param;

%% Load experiments 
load('/Users/ztakeo/Documents/GitHub/OED_ParamID/Codebase/SensAnalysis/max_sens_experiments_Tmax60.mat','results'); 

max_exp_num = results.max_exp_num_sorted;

% extract the unique experiment #'s while preserving the order ('stable' option)
max_exp_num_unique = unique(results.max_exp_num_sorted,'stable'); 
Num_unique_inputs = length(max_exp_num_unique);


input_path = 'InputLibrary/Experimental/Unformatted/';
%% Calculate Rc for each experiment
Rc = zeros(Num_unique_inputs,2);
Rc(:,1) = max_exp_num_unique;

selection_vector = [1;1;1;1;0;0;1;1;1;1;1;1;1;1;1;1;1;1;1;0;1;1;1;1;1]; %All, Year 2
SensFlag = 0; %Turn sensitivity calculation off or on
SensSelec = selection_vector;
sel_k = find(SensSelec);
Selected_params = theta_0(sel_k);
V_LM_CELL_sim = cell(Num_unique_inputs,1);
v_drop_norm = zeros(Num_unique_inputs,1);

for ii = 1:Num_unique_inputs
    % Reset R_c to nominal value
    run param/params_NCA
    
    % Load experimental data
    filename = strcat(input_path,num2str(max_exp_num_unique(ii)),'.mat');
    load(filename);
    
    % Simulate voltage response with nominal Rc value
    [V_LM_CELL_sim{ii},~] = DFN_sim_casadi(p, max_exp_num_unique(ii), Current_exp(1:11), Time_exp(1:11), V_LM_CELL(1:11), T_amb_sim(1:11), SensSelec, Selected_params, SensFlag); % SensFlag == 0

    % Plot before adjustment
%     figure
%     plot(Time_exp(1:11),V_LM_CELL(1:11))
%     hold on
%     plot(Time_exp(1:11),V_LM_CELL_sim{ii})
%     hold off
    
    % adjust Rc w/ Ohm's law and algebra
    I = Current_exp/p.Area; % DFN current convention
    Rc(ii,2) = p.R_c + ( (V_LM_CELL(11)-V_LM_CELL_sim{ii}(11)) - (V_LM_CELL(10)-V_LM_CELL_sim{ii}(10)) )/I(11);
    p.R_c = Rc(ii,2);
    
    % Resimulate 
    [V_LM_CELL_sim{ii},~] = DFN_sim_casadi(p, max_exp_num_unique(ii), Current_exp(1:11), Time_exp(1:11), V_LM_CELL(1:11), T_amb_sim(1:11), SensSelec, Selected_params, SensFlag); % SensFlag == 0
    v_drop_norm(ii) = norm(V_LM_CELL_sim{ii} - V_LM_CELL(1:11)); % calculate RMSE of drop for sanity check

    
    % Plot after adjustment
%     figure
%     plot(Time_exp(1:11),V_LM_CELL(1:11))
%     hold on
%     plot(Time_exp(1:11),V_LM_CELL_sim{ii})    
%     hold off
%     
%     disp('')
%     close all
end

