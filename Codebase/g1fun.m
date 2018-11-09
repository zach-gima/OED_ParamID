function [ output_args ] = g1fun( par0 )
%g1fun
%author: Dylan Kato
%Summary: takes parameters par0 a dim = 2 vector of the current iterations 
%g1 parameters and returns the scalar RMSE (for 5 g1 experiments)

%% Load parameter truth value and G1 experiments
verbose = false;
run param/params_truth
run param/params_NCA
% load g1 input
Inputs = load('m2m_comp/simulated_inputs/V_sim_G1.mat');

%% simulation set-up 
% Create parameter selection vector
SensSelec = [1;1;1;1;0;0;1;1;1;1;1;1;1;1;1;1;1;1;1;0;1]; %All params
sel_k = find(SensSelec);
G1ind = [3,4];
G2ind = [1,2,7,8,11,13];
% Set all parameter initial values to their true values
SelecParam = truth_param(sel_k); % Goes from 21x1 vector to 18x1
%THIS LINE IS KEY here we take the input parameters and put them into the
%vector that will be used for simulation
SelecParam(G1ind) = par0;

%%
% Parse experimental inputs
Current_exp = Inputs.Current_exp;
Time_exp = Inputs.Time_exp;
Voltage_exp = Inputs.V_LM_CELL;
y_dat = cell2mat(Voltage_exp);
SensFlag = 0; % no sensitivity calculations
num_inputs = length(Current_exp); % Determine # of inputs/experiments

%% This section was for my debugging
if verbose
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('#####################################')
    [SelecParam([3,4]), truth_param([3,4])]
    [SelecParam([3,4])./ truth_param([3,4])]
    disp('#####################################')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%

V = cell(num_inputs,1);
%parfor idx = 1:num_inputs
for idx = 1:num_inputs
    try
        V{idx} = DFN_sim_casadi(p,Current_exp{idx}, Time_exp{idx}, Voltage_exp{idx}, SensSelec, SelecParam, SensFlag);
        %right now I deal with min voltage erroy by truncating both the "truth"
        %and "simulated" voltages
        vv = Voltage_exp{idx};
        Voltage_exp{idx} = vv(1:length(V{idx}));
    catch j 
        if verbose
        disp('############')
        disp(['exp ', num2str(idx),' failed'])
        disp(j)
        disp('############')
        end
        %right now, if casadi fails, we just penalize heavily and arbitrarily
        %perhaps we could change this to a constraint to solve some issues!?!

        %set the simulated voltage for failed trials to zero
        V{idx} = zeros(length(Voltage_exp{idx}),1);
    end
end
y_dat = cell2mat(Voltage_exp);
%concatenate all of the voltage in a single column vector
V = cell2mat(V);
%calculate the output measurement (rmse)
output_args = sqrt(mean((y_dat - V).^2));
end