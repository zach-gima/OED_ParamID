%%% Script for tuning the IR drop by setting Rc. This is necessary for
%%% identifying parameters from experimental data
clearvars
close all
clc
fs = 25;

%% Load Parameters

run param/params_expID
theta_0 = expID_param;

%% Load experiments 
% load('/Users/ztakeo/Documents/GitHub/OED_ParamID/Codebase/SensAnalysis/max_sens_experiments_Tmax60.mat','results'); 
% 
% max_exp_num = results.max_exp_num_sorted;
% 
% % extract the unique experiment #'s while preserving the order ('stable' option)
% max_exp_num_unique = unique(results.max_exp_num_sorted,'stable'); 
% Num_unique_inputs = length(max_exp_num_unique);

% input_path = 'InputLibrary/Experimental/Unformatted/';
% input_path = 'InputLibrary/ValidationCycles/';
% input_path = '/Users/ztakeo/Desktop/OptExp/matfiles/';
input_path = '/Users/ztakeo/Documents/GitHub/OED_ParamID/Codebase/InputLibrary/ValidationCycles/Unformatted/';

% Num_unique_inputs = 12;
% max_exp_num_unique = [47; 57; 250; 270; 286; 287; 294; 298; 526; 632; 884; 898];

% Load sens files
val_files = dir(input_path);
% on Mac, use line below to ignore '.','..', and '.DS_Store' files that
% are loaded into r 
val_files = val_files(~ismember({val_files.name},{'.','..','.DS_Store'}));
Num_unique_inputs = length(val_files);

Rc_ID = cell(1,2);
for zz = 1:Num_unique_inputs
    max_exp_num_unique{zz,1} = val_files(zz).name;
    Rc_ID{zz,1} = max_exp_num_unique{zz};
end

%% Calculate Rc for each experiment
% Rc_ID = zeros(Num_unique_inputs,2);

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
    Rc_initial = p.R_c;

    % Load experimental data
%     filename = strcat(input_path,num2str(max_exp_num_unique(ii)),'.mat');
    filename = strcat(input_path,max_exp_num_unique{ii});
    load(filename);
    
%     V_LM_CELL = V;
    V_LM_CELL = Voltage_exp;
    
    %%% Supplement any missing variables
    if ~exist('Temp_exp','var')
        Temp_exp = ones(size(Voltage_exp))*25;
    elseif exist('Temp_exp','var')
        Temp_exp = Temp_exp{1};
    end
%     T_amb_sim = ones(Num_unique_inputs,1)*25;
    T_amb_sim = Temp_exp;
    
    % Add in Time
    if ~exist('Time_exp','var')
        Time_exp = (1:length(Voltage_exp))';
    end
    
    % Simulate voltage response with nominal Rc value
    [V_LM_CELL_sim{ii},~] = DFN_sim_casadi(p, max_exp_num_unique{ii}, Current_exp(1:30), Time_exp(1:30), V_LM_CELL(1:30), T_amb_sim(1:30), SensSelec, Selected_params, SensFlag,Rc_initial); % SensFlag == 0

%     Plot before adjustment
    figure('Position', [100 100 900 700])
    hold on

    plot(Time_exp(1:30),V_LM_CELL(1:30))
    plot(Time_exp(1:30),V_LM_CELL_sim{ii})
    legend('Experimental','Simulated')
    xlabel('Time (s)')
    ylabel('Voltage (V)')
    title('Before RC ID')
    set(gca,'Fontsize',fs);
    hold off
    
    % adjust Rc w/ Ohm's law and algebra 
    % Testing Profiles have 10s of 0 current, then the current profile
    % begins at t = 11s; therefore, we can calculate the IR drop by using
    % the voltage difference between t = 11 and t = 10s
    I = Current_exp/p.Area; % DFN current convention
    Rc_ID{ii,2} = p.R_c + ( (V_LM_CELL(11)-V_LM_CELL_sim{ii}(11)) - (V_LM_CELL(10)-V_LM_CELL_sim{ii}(10)) )/I(11);
    p.R_c = Rc_ID{ii,2};
    
    % Resimulate 
    [V_LM_CELL_sim{ii},~] = DFN_sim_casadi(p, max_exp_num_unique{ii}, Current_exp(1:30), Time_exp(1:30), V_LM_CELL(1:30), T_amb_sim(1:30), SensSelec, Selected_params, SensFlag,Rc_ID{ii,2}); % SensFlag == 0
    v_drop_norm(ii) = norm(V_LM_CELL_sim{ii} - V_LM_CELL(1:30)); % calculate RMSE of drop for sanity check
    
%     Plot after adjustment
    figure('Position', [100 100 900 700])
    hold on
    
    plot(Time_exp(1:30),V_LM_CELL(1:30))
    plot(Time_exp(1:30),V_LM_CELL_sim{ii})   
    legend('Experimental','Simulated')
    title('After RC ID')
    xlabel('Time (s)')
    ylabel('Voltage (V)')
    set(gca,'Fontsize',fs);
    hold off
    
    close all
    
    Rc = Rc_ID{ii,2};
    
    % Save Rc value to the experiment .mat file
    save(filename,'Current_exp','T_amb_sim','Time_exp','Voltage_exp','Rc')
    clear 'Current_exp' 'T_amb_sim' 'Time_exp' 'Voltage_exp' 'Rc' Temp_exp
end

