%% Script for pulling truth data from the sensitivity .mat files; don't need everything in the .mat file 
% and need the simulated data in a specific form
clearvars
close all
clc


%% User Input: File I/O 
% Directory location for sensitivity .mat files; *****make sure they only have
% the .mat files for the inputs in them.

% Uncomment Baseline to select

% Set path to sensitivity .mat files for every input in the library
% senspath = '/Users/ztakeo/Documents/GitHub/OED_ParamID/SensResults/Tmax60/';

%%%%% Baseline A: Full Parameter Set (1 Group)
% inputfinalpath = 'InputLibrary/MaxSensInputs/BaselineA_trim/';
% inputrawpath = strcat(inputfinalpath,'Unformatted/');
% load('/Users/ztakeo/Documents/GitHub/OED_ParamID/Codebase/SensAnalysis/max_sens_experiments_BaselineA_trim.mat','results'); 
% 
% % Set Number of Groups and params in each group
% Num_groups = 1; % desired number of param groups
% Group_size = 22; % for each group, specify # params to identify
% 
% % Create names for formatted inputs
% filename_input_vector{1} = strcat(inputfinalpath,'V_sim_G1.mat');

%%%%% Baseline B: Collinearity Only (1 Group)
% inputfinalpath = 'InputLibrary/MaxSensInputs/BaselineB_trim/';
% inputrawpath = strcat(inputfinalpath,'Unformatted/');
% load('/Users/ztakeo/Documents/GitHub/OED_ParamID/Codebase/SensAnalysis/max_sens_experiments_BaselineB_trim.mat','results'); 
% 
% % Set Number of Groups and params in each group
% Num_groups = 1; % desired number of param groups
% Group_size = 16; % for each group, specify # params to identify
% 
% % Create names for formatted inputs
% filename_input_vector{1} = strcat(inputfinalpath,'V_sim_G1.mat');

%%%% Baseline C: Collinearity + Sensitivity (2 Groups)
% inputfinalpath = 'InputLibrary/Experimental/';
inputfinalpath = 'InputLibrary/ValidationCycles/';
inputrawpath = strcat(inputfinalpath,'Unformatted/');
% load('/Users/ztakeo/Documents/GitHub/OED_ParamID/Codebase/SensAnalysis/max_sens_experiments_Tmax60_trim.mat','results'); 
% In experimental case, w/o results .mat file
% results.max_exp_num_sorted = [632;250;286;57;298;287;526;898;294;884;270;47;298];
% results.max_exp_num_sorted = [632;286;287;526;898;294;884;270;47];
% Load sens files
val_files = dir(inputrawpath);
% on Mac, use line below to ignore '.','..', and '.DS_Store' files that
% are loaded into r 
val_files = val_files(~ismember({val_files.name},{'.','..','.DS_Store'}));
Num_unique_inputs = length(val_files);

for zz = 1:Num_unique_inputs
    max_exp_num_unique{zz,1} = val_files(zz).name;
end
results.max_exp_num_sorted = max_exp_num_unique;
% Set Number of Groups and params in each group
Num_groups = 1; % desired number of param groups
Group_size = [24]; % for each group, specify # params to identify

% Create names for formatted inputs
filename_input_vector{1} = strcat(inputfinalpath,'validation_cycles.mat');
% filename_input_vector{1} = strcat(inputfinalpath,'V_sim_G1.mat');
% filename_input_vector{2} = strcat(inputfinalpath,'V_sim_G2.mat');
% filename_input_vector{3} = strcat(inputfinalpath,'V_sim_G3.mat');
% filename_input_vector{4} = strcat(inputfinalpath,'V_sim_G4.mat');
% filename_input_vector{5} = strcat(inputfinalpath,'V_sim_G2G1.mat');
% filename_input_vector{6} = strcat(inputfinalpath,'V_sim_G3G2G1.mat');
% filename_input_vector{7} = strcat(inputfinalpath,'V_sim_G4G3G2G1.mat');

%%%%% Perturbation Case
%%%%% Minus50
% inputfinalpath = 'InputLibrary/MaxSensInputs/minus50/';
% inputrawpath = strcat(inputfinalpath,'Unformatted/');
% load('/Users/ztakeo/Documents/GitHub/OED_ParamID/Codebase/SensAnalysis/max_sens_experiments_minus50.mat','results'); 
% 
% % Set Number of Groups and params in each group
% Num_groups = 2; % desired number of param groups
% Group_size = [6,7]; % for each group, specify # params to identify
% 
% % Create names for formatted inputs
% filename_input_vector{1} = strcat(inputfinalpath,'V_sim_G1.mat');
% filename_input_vector{2} = strcat(inputfinalpath,'V_sim_G2.mat');
% filename_input_vector{3} = strcat(inputfinalpath,'V_sim_G2G1.mat');

%%%%% Plus50
% inputfinalpath = 'InputLibrary/MaxSensInputs/plus50/';
% inputrawpath = strcat(inputfinalpath,'Unformatted/');
% load('/Users/ztakeo/Documents/GitHub/OED_ParamID/Codebase/SensAnalysis/max_sens_experiments_plus50.mat','results'); 
% 
% % Set Number of Groups and params in each group
% Num_groups = 2; % desired number of param groups
% Group_size = [6,7]; % for each group, specify # params to identify
% 
% % Create names for formatted inputs
% filename_input_vector{1} = strcat(inputfinalpath,'V_sim_G1.mat');
% filename_input_vector{2} = strcat(inputfinalpath,'V_sim_G2.mat');
% filename_input_vector{3} = strcat(inputfinalpath,'V_sim_G2G1.mat');

% if exist(inputrawpath,'dir') == 7
%     rmdir(inputrawpath,'s'); % delete folder first (assuming it exists); this prevents the folder from keeping older .mat files
% end
% mkdir(inputrawpath);

%% Load #s of the Max Sens Inputs selected

% Pull out the #'s for the max sens. experiments
max_exp_num = results.max_exp_num_sorted;
Num_params = length(max_exp_num);
Num_inputs = Num_params;

% extract the unique experiment #'s while preserving the order ('stable' option)
max_exp_num_unique = unique(results.max_exp_num_sorted,'stable'); 
Num_unique_inputs = length(max_exp_num_unique);

% % Copy Max Sens Inputs into new directory
% for ii = 1:Num_unique_inputs
%     filename = strcat(senspath,num2str(max_exp_num_unique(ii)),'.mat');
%     copyfile(filename,inputrawpath)
% end

%% Quick Error Checking
if length(Group_size) ~= Num_groups
    error('# groups specified in Group_size does not match Group_size. Fix discrepancy')
end
% if sum(Group_size) ~= Num_params
%     error('Total # of params specified in Group_size does not match Num_params to be identified.')
% end

%% Concatenate Input and Corresponding Simulation Data Across Group Inputs

% exp_num_track = [];
exp_idx = 1; % index variable used to iterate through all max_exp_num entries

for mm = 1:Num_groups
    
    % Initialize variables
   
    Current_exp_cell = {};
    Time_exp_cell = {};
    V_LM_CELL = {};
    cssn_sim_cell = {};
    cssp_sim_cell = {};
    etan_sim_cell = {};
    etap_sim_cell = {};
    T1_sim_cell = {};
    T2_sim_cell = {};
    T_amb_sim_cell = {};
    exp_num_cell = {};
    Rc_exp_cell = {};
    
    for zz = 1:Group_size(mm) 
        % Check whether input has already been formatted; sometimes same input
        % is selected by multiple parameters; we only want to run once
%         if ismember(max_exp_num(zz),exp_num_track)
%             continue;
%         else
%         end
%         % Put exp_num into the exp num tracking array so that it won't be
%         % duplicated
%         exp_num_track = vertcat(exp_num_track,max_exp_num(exp_idx));

        % Load experiment
%         exp_str = strcat(num2str(max_exp_num(exp_idx)),'.mat');
        exp_str = max_exp_num{exp_idx};

%         load(strcat(inputrawpath,exp_str),'Current_exp','Time_exp','V','alg_states','T_amb');
        load(strcat(inputrawpath,exp_str));

        % I, V, t
        Current_exp_cell = vertcat(Current_exp_cell,Current_exp);
        Time_exp_cell = vertcat(Time_exp_cell,Time_exp);
        V_LM_CELL = vertcat(V_LM_CELL,Voltage_exp);
        Rc_exp_cell = vertcat(Rc_exp_cell,Rc);

        % algebraic states
%         cssn_sim_cell = vertcat(cssn_sim_cell,alg_states.cssn_sim');
%         cssp_sim_cell = vertcat(cssp_sim_cell,alg_states.cssp_sim');
%         etan_sim_cell = vertcat(etan_sim_cell,alg_states.etan_sim');
%         etap_sim_cell = vertcat(etap_sim_cell,alg_states.etap_sim');

        % temperature
%         T1_sim_cell = vertcat(T1_sim_cell,alg_states.T1_sim');
%         T2_sim_cell = vertcat(T2_sim_cell,alg_states.T2_sim');
        T_amb_sim_cell = vertcat(T_amb_sim_cell,T_amb_sim(1));
        
        % exp number
        exp_num_cell = vertcat(exp_num_cell,num2str(max_exp_num{exp_idx}));
        
        % Min V checking (issue that comes up during model-to-model comp)
        min_V = min(V_LM_CELL{end});
        fprintf('Exp. %s has min V = %0.3f \n',max_exp_num{exp_idx},min_V);
        
        max_V = max(V_LM_CELL{end});
        fprintf('Exp. %s has max V = %0.3f \n',max_exp_num{exp_idx},max_V);
        % Update exp_idx, used to keep track of which max_exp_num we're on
        exp_idx = exp_idx + 1;
    end
    
    % rename the variables and remove _cell (just so it plays nicely w/ Param_ID.m and DFN_sim_casadi.m) 
    Current_exp = Current_exp_cell;
    Time_exp = Time_exp_cell;
    T_amb_sim = T_amb_sim_cell;
    exp_num = exp_num_cell;
    Rc = Rc_exp_cell;
    
    % If including algebraic states
%     cssn_sim = cssn_sim_cell;
%     cssp_sim = cssp_sim_cell;
%     etan_sim = etan_sim_cell;
%     etap_sim = etap_sim_cell;
%     T1_sim = T1_sim_cell;
%     T2_sim = T2_sim_cell;
    
%     save(filename_input_vector{mm},'Time_exp','Current_exp','V_LM_CELL','cssn_sim','cssp_sim','etan_sim','etap_sim','T1_sim','T2_sim','T_amb_sim','exp_num');
    save(filename_input_vector{mm},'Time_exp','Current_exp','V_LM_CELL','T_amb_sim','exp_num','Rc');

end

%% Uncomment to create .mat files where groups are combined G2G1, G3G2G1...
%%% Pulled from DFN_sim_main
% 
% %%% Group 1 Input Data
% S1 = load(filename_input_vector{1});
% Current_exp_1 = S1.Current_exp;
% Time_exp_1 = S1.Time_exp;
% V_LM_CELL_1 = S1.V_LM_CELL;
% T_amb_sim_1 = S1.T_amb_sim;
% exp_num_1 = S1.exp_num;
% Rc_1 = S1.Rc;
% 
% % % If including algebraic states
% % cssn_sim_1 = S1.cssn_sim;
% % cssp_sim_1 = S1.cssp_sim;
% % etan_sim_1 = S1.etan_sim;
% % etap_sim_1 = S1.etap_sim;
% % T1_sim_1 = S1.T1_sim;
% % T2_sim_1 = S1.T2_sim;
% 
% %%% Group 2 Input Data
% S2 = load(filename_input_vector{2});
% Current_exp_2 = S2.Current_exp;
% Time_exp_2 = S2.Time_exp;
% V_LM_CELL_2 = S2.V_LM_CELL;
% T_amb_sim_2 = S2.T_amb_sim;
% exp_num_2 = S2.exp_num;
% Rc_2 = S2.Rc;
% 
% % % If including algebraic states
% % cssn_sim_2 = S2.cssn_sim;
% % cssp_sim_2 = S2.cssp_sim;
% % etan_sim_2 = S2.etan_sim;
% % etap_sim_2 = S2.etap_sim;
% % T1_sim_2 = S2.T1_sim;
% % T2_sim_2 = S2.T2_sim;
% 
% 
% %%% Group 3 Input Data
% S3 = load(filename_input_vector{3});
% Current_exp_3 = S3.Current_exp;
% Time_exp_3 = S3.Time_exp;
% V_LM_CELL_3 = S3.V_LM_CELL;
% T_amb_sim_3 = S3.T_amb_sim;
% exp_num_3 = S3.exp_num;
% Rc_3 = S3.Rc;
% 
% 
% %%% Group 4 Input Data
% S4 = load(filename_input_vector{4});
% Current_exp_4 = S4.Current_exp;
% Time_exp_4 = S4.Time_exp;
% V_LM_CELL_4 = S4.V_LM_CELL;
% T_amb_sim_4 = S4.T_amb_sim;
% exp_num_4 = S4.exp_num;
% Rc_4 = S4.Rc;
% 
% %%% G2G1
% Current_exp =  vertcat(Current_exp_2,Current_exp_1);
% Time_exp =  vertcat(Time_exp_2,Time_exp_1);
% V_LM_CELL = vertcat(V_LM_CELL_2,V_LM_CELL_1);
% T_amb_sim = vertcat(T_amb_sim_2,T_amb_sim_1);
% exp_num = vertcat(exp_num_2,exp_num_1);
% Rc = vertcat(Rc_2,Rc_1);
% 
% % If including algebraic states
% % cssn_sim = vertcat(cssn_sim_2,cssn_sim_1);
% % cssp_sim = vertcat(cssp_sim_2,cssp_sim_1);
% % etan_sim = vertcat(etan_sim_2,etan_sim_1);
% % etap_sim = vertcat(etap_sim_2,etap_sim_1);
% % T1_sim = vertcat(T1_sim_2,T1_sim_1);
% % T2_sim = vertcat(T2_sim_2,T2_sim_1);
% 
% save(filename_input_vector{5},'Time_exp','Current_exp','V_LM_CELL','T_amb_sim','exp_num','Rc');
% 
% %%% G3G2G1
% Current_exp =  vertcat(Current_exp_3,Current_exp_2,Current_exp_1);
% Time_exp =  vertcat(Time_exp_3,Time_exp_2,Time_exp_1);
% V_LM_CELL = vertcat(V_LM_CELL_3,V_LM_CELL_2,V_LM_CELL_1);
% T_amb_sim = vertcat(T_amb_sim_3,T_amb_sim_2,T_amb_sim_1);
% exp_num = vertcat(exp_num_3,exp_num_2,exp_num_1);
% Rc = vertcat(Rc_3,Rc_2,Rc_1);
% 
% save(filename_input_vector{6},'Time_exp','Current_exp','V_LM_CELL','T_amb_sim','exp_num','Rc');
% 
% 
% %%% G4G3G2G1
% Current_exp =  vertcat(Current_exp_4,Current_exp_3,Current_exp_2,Current_exp_1);
% Time_exp =  vertcat(Time_exp_4,Time_exp_3,Time_exp_2,Time_exp_1);
% V_LM_CELL = vertcat(V_LM_CELL_4,V_LM_CELL_3,V_LM_CELL_2,V_LM_CELL_1);
% T_amb_sim = vertcat(T_amb_sim_4,T_amb_sim_3,T_amb_sim_2,T_amb_sim_1);
% exp_num = vertcat(exp_num_4,exp_num_3,exp_num_2,exp_num_1);
% Rc = vertcat(Rc_4,Rc_3,Rc_2,Rc_1);
% 
% save(filename_input_vector{7},'Time_exp','Current_exp','V_LM_CELL','T_amb_sim','exp_num','Rc');
% 
% % save(filename_input_vector{3},'Current_exp','Time_exp','V_LM_CELL','cssn_sim','cssp_sim','etan_sim','etap_sim','T1_sim','T2_sim','T_amb_sim','exp_num')
% % save(filename_input_vector{3},'Time_exp','Current_exp','V_LM_CELL','T_amb_sim','exp_num','Rc');
