%% Confidence Interval After the Fact
% Need to re-calculate confidence intervals after ParamID finished
clear all
close all
clc

%% User Inputs
run param/params_NCA
num_groups = 4; %4 for nominal case

results_folder = 'm2m_comp/ID_results/02-Aug-2018 07_57_02/';
input_folder = 'm2m_comp/simulated_inputs/';

%% Load Variables
% Setup ID'ed param cell 
park0_input_cell = cell(4,1);

park0_input_cell{1} = strcat(results_folder,'G1_nom.mat');
park0_input_cell{2} = strcat(results_folder,'G2G1_nom.mat');
park0_input_cell{3} = strcat(results_folder,'G3G2G1_nom.mat');
park0_input_cell{4} = strcat(results_folder,'G4G3G2G1_nom.mat');

selection_vector = zeros(21,4); %21 parameters
selection_vector(:,1) = [0;0;1;1;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0]; %G1
selection_vector(:,2) = [1;1;1;1;0;0;0;0;1;1;0;0;1;0;1;0;0;0;0;0;0]; %G2
selection_vector(:,3) = [1;1;1;1;0;0;0;0;1;1;0;1;1;0;1;1;0;1;1;0;1]; %G3
selection_vector(:,4) = [1;1;1;1;0;0;1;1;1;1;1;1;1;1;1;1;1;1;1;0;1]; %G4

% Create vector of indices for calculating confidence intervals
% Will be calculating confidence intervals separately for each group
ci_select = cell(4,1);
ci_select{1} = find(selection_vector(:,1));
ci_select{2} = find(selection_vector(:,2) - selection_vector(:,1));
ci_select{3} = find(selection_vector(:,3) - selection_vector(:,2));
ci_select{4} = find(selection_vector(:,4) - selection_vector(:,3));

% CI Input filenames
% Note: for calculating confidence intervals, will only use the oed-cvx
% selected inputs for that respective group
ci_input_vector = cell(4,1);
ci_input_vector{1} = strcat(input_folder,'V_sim_G1.mat');
ci_input_vector{2} = strcat(input_folder,'V_sim_G2.mat');
ci_input_vector{3} = strcat(input_folder,'V_sim_G3.mat');
ci_input_vector{4} = strcat(input_folder,'V_sim_G4.mat');

ci95_full = zeros(21,1); %vector for storing confidence interval for params

for ii = 1:num_groups
    % load results from ParamID -- parse out identified param values
    ID_results = load(park0_input_cell{ii}); % loads a bunch of variables; importantly: paramID_out
    paramID_out = ID_results.paramID_out;
    
    %load inputs
    ci_inputs = load(ci_input_vector{ii});
    
    % Extract identified params from data file
    park0 = paramID_out.save_param_org(:,end); % grad final set of identified params
    
    SensSelec = selection_vector(:,ii);

    % Calculate local sensitivity and final RMSE
    [J_LM,ci95,sigma_y,covar_p] = conf_interval(p, ci_inputs, SensSelec, park0);

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