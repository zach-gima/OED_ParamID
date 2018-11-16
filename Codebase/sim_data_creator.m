%% Script for pulling truth data from the sensitivity .mat files; don't need everything in the .mat file 
% and need the simulated data in a specific form
clearvars
close all
clc

%% Load #s of the Max Sens Inputs selected

%% Copy

%% Specify Input File Directory
input_path = '/Users/ztakeo/Documents/GitHub/DFN_CasADi/3_Journal_JPS/InputLibrary/MaxSensInputs/G1/';

%% G1: Load All Inputs for Group 1 Parameters (Highest sens. G1 param to
%%% lowest)
Num_inputs = 7;

load(strcat(input_path,'610.mat'),'Current_exp','Time_exp','V','alg_states','T_amb');

% I, V, t
Current_exp1 = Current_exp';
Time_exp1 = Time_exp';
V_LM_CELL1 = V;
% Rc_exp1 = Rc;

% algebraic states
cssn_sim1 = alg_states.cssn_sim';
cssp_sim1 = alg_states.cssp_sim';
etan_sim1 = alg_states.etan_sim';
etap_sim1 = alg_states.etap_sim';

% temperature
T1_sim1 = alg_states.T1_sim';
T2_sim1 = alg_states.T2_sim';
T_amb_sim1 = T_amb(1);

load(strcat(input_path,'83.mat'),'Current_exp','Time_exp','V','alg_states','T_amb');

% I, V, t
Current_exp2 = Current_exp';
Time_exp2 = Time_exp';
V_LM_CELL2 = V;
% Rc_exp1 = Rc;

% algebraic states
cssn_sim2 = alg_states.cssn_sim';
cssp_sim2 = alg_states.cssp_sim';
etan_sim2 = alg_states.etan_sim';
etap_sim2 = alg_states.etap_sim';

% temperature
T1_sim2 = alg_states.T1_sim';
T2_sim2 = alg_states.T2_sim';
T_amb_sim2 = T_amb(1);

load(strcat(input_path,'867.mat'),'Current_exp','Time_exp','V','alg_states','T_amb');

% I, V, t
Current_exp3 = Current_exp';
Time_exp3 = Time_exp';
V_LM_CELL3 = V;
% Rc_exp1 = Rc;

% algebraic states
cssn_sim3 = alg_states.cssn_sim';
cssp_sim3 = alg_states.cssp_sim';
etan_sim3 = alg_states.etan_sim';
etap_sim3 = alg_states.etap_sim';

% temperature
T1_sim3 = alg_states.T1_sim';
T2_sim3 = alg_states.T2_sim';
T_amb_sim3 = T_amb(1);

load(strcat(input_path,'875.mat'),'Current_exp','Time_exp','V','alg_states','T_amb');

% I, V, t
Current_exp4 = Current_exp';
Time_exp4 = Time_exp';
V_LM_CELL4 = V;
% Rc_exp1 = Rc;

% algebraic states
cssn_sim4 = alg_states.cssn_sim';
cssp_sim4 = alg_states.cssp_sim';
etan_sim4 = alg_states.etan_sim';
etap_sim4 = alg_states.etap_sim';

% temperature
T1_sim4 = alg_states.T1_sim';
T2_sim4 = alg_states.T2_sim';
T_amb_sim4 = T_amb(1);

load(strcat(input_path,'202.mat'),'Current_exp','Time_exp','V','alg_states','T_amb');

% I, V, t
Current_exp5 = Current_exp';
Time_exp5 = Time_exp';
V_LM_CELL5 = V;
% Rc_exp1 = Rc;

% algebraic states
cssn_sim5 = alg_states.cssn_sim';
cssp_sim5 = alg_states.cssp_sim';
etan_sim5 = alg_states.etan_sim';
etap_sim5 = alg_states.etap_sim';

% temperature
T1_sim5 = alg_states.T1_sim';
T2_sim5 = alg_states.T2_sim';
T_amb_sim5 = T_amb(1);

load(strcat(input_path,'567.mat'),'Current_exp','Time_exp','V','alg_states','T_amb');

% I, V, t
Current_exp6 = Current_exp';
Time_exp6 = Time_exp';
V_LM_CELL6 = V;
% Rc_exp1 = Rc;

% algebraic states
cssn_sim6 = alg_states.cssn_sim';
cssp_sim6 = alg_states.cssp_sim';
etan_sim6 = alg_states.etan_sim';
etap_sim6 = alg_states.etap_sim';

% temperature
T1_sim6 = alg_states.T1_sim';
T2_sim6 = alg_states.T2_sim';
T_amb_sim6 = T_amb(1);

load(strcat(input_path,'1247.mat'),'Current_exp','Time_exp','V','alg_states','T_amb');

% I, V, t
Current_exp7 = Current_exp';
Time_exp7 = Time_exp';
V_LM_CELL7 = V;
% Rc_exp1 = Rc;

% algebraic states
cssn_sim7 = alg_states.cssn_sim';
cssp_sim7 = alg_states.cssp_sim';
etan_sim7 = alg_states.etan_sim';
etap_sim7 = alg_states.etap_sim';

% temperature
T1_sim7 = alg_states.T1_sim';
T2_sim7 = alg_states.T2_sim';
T_amb_sim7 = T_amb(1);


%% G2: Load All Inputs for Group 2 Parameters (Highest sens. G1 param to
%%% lowest)

Finish for G2 

%% Concatenate Data Across All Inputs in Each Group
%%% Since truth data already generated during sensitivity calculations,
%%% just need to format it properly

% Each variable should be a cell (size = # inputs x 1)
Current_exp = cell(Num_inputs,1);
Time_exp = cell(Num_inputs,1);
V_LM_CELL = cell(Num_inputs,1);
% Rc_exp = cell(Num_inputs,1);

cssn_sim = cell(Num_inputs,1);
cssp_sim = cell(Num_inputs,1);
etan_sim = cell(Num_inputs,1);
etap_sim = cell(Num_inputs,1);

T1_sim = cell(Num_inputs,1);
T2_sim = cell(Num_inputs,1);
T_amb_sim = cell(Num_inputs,1);

% Concatenate Input and Corresponding Simulation Data Across Group Inputs
for i=1:Num_inputs 
    Current_exp{i} = eval(['Current_exp' num2str(i)]);
    Time_exp{i} = eval(['Time_exp' num2str(i)]);
    V_LM_CELL{i} = eval(['V_LM_CELL' num2str(i)]);
%     Rc_exp{i} = eval(['Rc_exp' num2str(i)]);

    cssn_sim{i} = eval(['cssn_sim' num2str(i)]);
    cssp_sim{i} = eval(['cssp_sim' num2str(i)]);
    etan_sim{i} = eval(['etan_sim' num2str(i)]);
    etap_sim{i} = eval(['etap_sim' num2str(i)]);

    T1_sim{i} = eval(['T1_sim' num2str(i)]);
    T2_sim{i} = eval(['T2_sim' num2str(i)]);
    T_amb_sim{i} = eval(['T_amb_sim' num2str(i)]);

end

save('V_sim_G1.mat','Time_exp','Current_exp','V_LM_CELL','cssn_sim','cssp_sim','etan_sim','etap_sim','T1_sim','T2_sim','T_amb_sim');
%save('V_sim_G2.mat','Time_exp','Current_exp','V_LM_CELL','cssn_sim','cssp_sim','etan_sim','etap_sim','T1_sim','T2_sim','T_amb_sim');


%% Uncomment to create .mat files where groups are combined G2G1, G3G2G1...
%%% Pulled from DFN_sim_main

% S1 = load('V_sim_G1.mat');
% Current_exp_1 = S1.Current_exp;
% Time_exp_1 = S1.Time_exp;
% V_LM_CELL_1 = S1.V_LM_CELL;
% cssn_sim_1 = S1.cssn_sim;
% cssp_sim_1 = S1.cssp_sim;
% etan_sim_1 = S1.etan_sim;
% etap_sim_1 = S1.etap_sim;
% T1_sim_1 = S1.T1_sim;
% T2_sim_1 = S1.T2_sim;
% T_amb_sim_1 = S1.T_amb_sim;
% 
% S2 = load('V_sim_G2.mat');
% Current_exp_2 = S2.Current_exp;
% Time_exp_2 = S2.Time_exp;
% V_LM_CELL_2 = S2.V_LM_CELL;
% cssn_sim_2 = S2.cssn_sim;
% cssp_sim_2 = S2.cssp_sim;
% etan_sim_2 = S2.etan_sim;
% etap_sim_2 = S2.etap_sim;
% T1_sim_2 = S1.T1_sim;
% T2_sim_2 = S1.T2_sim;
% T_amb_sim_2 = S1.T_amb_sim;

% % G2G1
% Current_exp =  vertcat(Current_exp_2,Current_exp_1);
% Time_exp =  vertcat(Time_exp_2,Time_exp_1);
% V_LM_CELL = vertcat(V_LM_CELL_2,V_LM_CELL_1);
% cssn_sim = vertcat(cssn_sim_2,cssn_sim_1);
% cssp_sim = vertcat(cssp_sim_2,cssp_sim_1);
% etan_sim = vertcat(etan_sim_2,etan_sim_1);
% etap_sim = vertcat(etap_sim_2,etap_sim_1);
% T1_sim = vertcat(T1_sim_2,T1_sim_1);
% T2_sim = vertcat(T2_sim_2,T2_sim_1);
% T_amb_sim = vertcat(T_amb_sim_2,T_amb_sim_1);
% 
% save('V_sim_G2G1.mat','Current_exp','Time_exp','V_LM_CELL','cssn_sim','cssp_sim','etan_sim','etap_sim')