%% Confidence Interval After the Fact
% Need to re-calculate confidence intervals after ParamID finished
clear all
close all
clc

%% User Inputs
run param/params_NCA

% ParamID Initial Conditions
run param/params_nominal
run param/params_bounds
theta_0 = Nominal_param;
% run param/params_truth % loads truth_param array
num_groups = 1; %4 for nominal case

output_folder = 'C:/Users/zgima/Box Sync/HPC/HPC2/B+5/';
park0_input_cell{1} = strcat(output_folder,'collinearity_nom.mat');
filename_output = strcat(output_folder,'collinearity_nom.mat');

input_folder = 'C:/Users/zgima/Documents/GitHub/OED_ParamID/Codebase/InputLibrary/MaxSensInputs/BaselineB_trim/';
% CI Input filenames
% Note: for calculating confidence intervals, will only use the oed-cvx
% selected inputs for that respective group
ci_input_vector{1} = strcat(input_folder,'V_sim_G1.mat');

% t_paramID = [0;72*3600];

% Set Levenberg-Marquardt Conditions
param_exit_thresh = 1e-3; % 1e-6 default values used in matlab (called StepTol)
chi_sq_exit_thresh = 1e-6; % 1e-6 default values used in matlab (called FuncTol)
chi_sq_Abs_exit_thresh = 1e-4; % Absolute cost function tolerance (Func Absolute Tol)

LM_options.exit_cond = [param_exit_thresh, chi_sq_exit_thresh, chi_sq_Abs_exit_thresh];
LM_options.maxIter = 20;

%% Load Variables
% Setup ID'ed param cell 
% selection_vector(:,1) =
% [1;1;1;1;0;0;1;1;1;1;1;1;1;1;1;1;1;1;1;0;1;1;1;1;1]; % Baseline A
selection_vector(:,1) = [1;1;1;1;0;0;1;1;1;0;1;1;0;1;1;0;0;1;0;0;1;1;1;1;0]; % Baseline B

%%% For perturbation analysis
% For perturbation analysis specifically where perturbing away from nominal
truth_param = Nominal_param;

% Perturb parameters of interest all at the beginning
perturb_index = find(selection_vector(:,1)); % G1
% perturb_index = find(selection_vector(:,2)); %G1 & G2

perturb_factor = 1.05;
theta_0(perturb_index) = perturb_factor*theta_0(perturb_index);
theta_0_true = theta_0;

% Create vector of indices for calculating confidence intervals
% Will be calculating confidence intervals separately for each group
ci_select{1} = find(selection_vector(:,1));
ci95_full = zeros(25,1); %vector for storing confidence interval for params
rmse_final = [];

for ii = 1:num_groups
    % load results from ParamID -- parse out identified param values
    load(park0_input_cell{ii}); % loads a bunch of variables; importantly: paramID_out
    paramID_out = ID_results.paramID_out;
    
    iter_history(ii) = length(paramID_out.save_param_org);
    %load inputs
    ci_inputs = load(ci_input_vector{ii});
    
    % Extract identified params from data file
    park0 = paramID_out.save_param_org(:,end); % grad final set of identified params
    
    SensSelec = selection_vector(:,ii);

    tic
    % Calculate local sensitivity and final RMSE
    [J_LM,ci95,sigma_y,covar_p,alg_states] = conf_interval(p, ci_inputs, SensSelec, park0);
    
    t_ci = toc;
    
    % Save very initial RMSE
    if isempty(rmse_final)
        rmse_final = vertcat(rmse_final,paramID_out.save_RMSE(1));
    end
    
    %Append RMSE for current group
    rmse_final = vertcat(rmse_final,paramID_out.save_RMSE(end));

    % Store Confidence Intervals
    sel_k = find(SensSelec);
    ci95 = horzcat(ci95,sel_k); % append indices for each parameter; will be used below
    
    % After the 1st group, need to only replace the confidence intervals
    % calculated for the current group. This is necessary because
    % conf_interval calculates for all the identified parameters passed in.
    % For example, in G2G1 group, G1 and G2 conf intervals are calcualted,
    % but we've already calculated the correct G1 intervals the previous
    % iteration and don't want to overwrite them.
    if ii > 1
        %iterate through the indices in ci_select
        for kk = 1:length(ci_select{ii})
            % iterate through the indices in sel_k, the 2nd column of ci95
            for mm = 1:size(ci95,1)
                % when the sel_k index matches the ci_select index, then
                % write that conf. interval value to ci95_full, which
                % stores everything
                if ci95(mm,2) == ci_select{ii}(kk)
                    ci95_full(ci_select{ii}(kk)) = ci95(mm,1); 
                    break;
                end
            end
        end
    else
        ci95_full(ci_select{1}) = ci95(:,1);
    end
end

%%%% SAVE THE CONFIDENCE INTERVAL BEFORE MESSING WITH THIS
save(filename_output,'paramID_out','rmse_final','ci95_full','sigma_y',...
'covar_p','t_paramID','t_ci',...
'iter_history','output_folder','sel_k','selection_vector','truth_param',...
'theta_0_true','LM_options','J_LM','bounds','alg_states');

Param_ID_plot(truth_param,theta_0_true,sel_k,paramID_out,ci95_full,t_paramID,rmse_final,output_folder,LM_options,bounds,alg_states,selection_vector)

disp('Confidence interval calculation complete.')