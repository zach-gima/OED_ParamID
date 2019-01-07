inputrawpath = 'InputLibrary/ValidationCycles/Unformatted/';

%% Concatenate Input and Corresponding Simulation Data Across Group Inputs
Num_groups = 1;
% exp_num_track = [];
exp_idx = 1; % index variable used to iterate through all val_exp_num entries

val_exp_num = ['C1';'C2';'C3';'C4';'C5';'C6';'C7';'C8';'D1';'D2';'D3';'D4';'D5';'D6';'D7';'D8'];

for mm = 1:Num_groups
    
    % Initialize variables
   
    Current_exp_cell = {};
    Time_exp_cell = {};
    V_LM_CELL = {};
    T_amb_sim_cell = {};
    exp_num_cell = {};
    for zz = 1:Group_size(mm) 
        % Check whether input has already been formatted; sometimes same input
        % is selected by multiple parameters; we only want to run once
%         if ismember(val_exp_num(zz),exp_num_track)
%             continue;
%         else
%         end
%         % Put exp_num into the exp num tracking array so that it won't be
%         % duplicated
%         exp_num_track = vertcat(exp_num_track,val_exp_num(exp_idx));

        % Load experiment
        exp_str = strcat(val_exp_num(exp_idx),'.mat');
        load(strcat(inputrawpath,exp_str),'Current_exp','Time_exp','V','alg_states','T_amb');

        % I, V, t
        Current_exp_cell = vertcat(Current_exp_cell,Current_exp');
        Time_exp_cell = vertcat(Time_exp_cell,Time_exp');
        V_LM_CELL = vertcat(V_LM_CELL,V);
        % Rc_exp_cell = Rc;

        % algebraic states
        cssn_sim_cell = vertcat(cssn_sim_cell,alg_states.cssn_sim');
        cssp_sim_cell = vertcat(cssp_sim_cell,alg_states.cssp_sim');
        etan_sim_cell = vertcat(etan_sim_cell,alg_states.etan_sim');
        etap_sim_cell = vertcat(etap_sim_cell,alg_states.etap_sim');

        % temperature
        T1_sim_cell = vertcat(T1_sim_cell,alg_states.T1_sim');
        T2_sim_cell = vertcat(T2_sim_cell,alg_states.T2_sim');
        T_amb_sim_cell = vertcat(T_amb_sim_cell,T_amb(1));
        
        % exp number
        exp_num_cell = vertcat(exp_num_cell,num2str(val_exp_num(exp_idx)));
        
        % Debugging
        if length(Time_exp') ~= length(alg_states.T1_sim')
            fprintf('Experiment %s has mismatched data sizes \n',exp_str);
        end
        
        % Min V checking (issue that comes up during model-to-model comp)
        min_V = min(V_LM_CELL{end});
        fprintf('Exp. %i has min V = %0.3f \n',val_exp_num(exp_idx),min_V);
        
        max_V = max(V_LM_CELL{end});
        fprintf('Exp. %i has max V = %0.3f \n',val_exp_num(exp_idx),max_V);
        % Update exp_idx, used to keep track of which val_exp_num we're on
        exp_idx = exp_idx + 1;
    end
    
    % rename the variables and remove _cell (just so it plays nicely w/ Param_ID.m and DFN_sim_casadi.m) 
    Current_exp = Current_exp_cell;
    Time_exp = Time_exp_cell;
    cssn_sim = cssn_sim_cell;
    cssp_sim = cssp_sim_cell;
    etan_sim = etan_sim_cell;
    etap_sim = etap_sim_cell;
    T1_sim = T1_sim_cell;
    T2_sim = T2_sim_cell;
    T_amb_sim = T_amb_sim_cell;
    exp_num = exp_num_cell;
    
    save(filename_input_vector{mm},'Time_exp','Current_exp','V_LM_CELL','cssn_sim','cssp_sim','etan_sim','etap_sim','T1_sim','T2_sim','T_amb_sim','exp_num');
end

%% Uncomment to create .mat files where groups are combined G2G1, G3G2G1...
%%% Pulled from DFN_sim_main

S1 = load(filename_input_vector{1});
Current_exp_1 = S1.Current_exp;
Time_exp_1 = S1.Time_exp;
V_LM_CELL_1 = S1.V_LM_CELL;
cssn_sim_1 = S1.cssn_sim;
cssp_sim_1 = S1.cssp_sim;
etan_sim_1 = S1.etan_sim;
etap_sim_1 = S1.etap_sim;
T1_sim_1 = S1.T1_sim;
T2_sim_1 = S1.T2_sim;
T_amb_sim_1 = S1.T_amb_sim;
exp_num_1 = S1.exp_num;

S2 = load(filename_input_vector{2});
Current_exp_2 = S2.Current_exp;
Time_exp_2 = S2.Time_exp;
V_LM_CELL_2 = S2.V_LM_CELL;
cssn_sim_2 = S2.cssn_sim;
cssp_sim_2 = S2.cssp_sim;
etan_sim_2 = S2.etan_sim;
etap_sim_2 = S2.etap_sim;
T1_sim_2 = S2.T1_sim;
T2_sim_2 = S2.T2_sim;
T_amb_sim_2 = S2.T_amb_sim;
exp_num_2 = S2.exp_num;


% G2G1
Current_exp =  vertcat(Current_exp_2,Current_exp_1);
Time_exp =  vertcat(Time_exp_2,Time_exp_1);
V_LM_CELL = vertcat(V_LM_CELL_2,V_LM_CELL_1);
cssn_sim = vertcat(cssn_sim_2,cssn_sim_1);
cssp_sim = vertcat(cssp_sim_2,cssp_sim_1);
etan_sim = vertcat(etan_sim_2,etan_sim_1);
etap_sim = vertcat(etap_sim_2,etap_sim_1);
T1_sim = vertcat(T1_sim_2,T1_sim_1);
T2_sim = vertcat(T2_sim_2,T2_sim_1);
T_amb_sim = vertcat(T_amb_sim_2,T_amb_sim_1);
exp_num = vertcat(exp_num_2,exp_num_1);
save(filename_input_vector{3},'Current_exp','Time_exp','V_LM_CELL','cssn_sim','cssp_sim','etan_sim','etap_sim','T1_sim','T2_sim','T_amb_sim','exp_num')