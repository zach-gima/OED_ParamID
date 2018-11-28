    %% Created by Dylan Kato, Zach Gima 2018-11-15
% This script first reads in pre-computed sensitivities for a pre-defined set of parameters in an electrochemical battery model. 
% It then clusters and eliminates parameters by analyzing the collinearity
% among the parameters and sensitivity magnitude of each parameter across
% all inputs in the library. Finally, the script determines the "optimal" inputs to use for identifying the reduced set of parameters 
% -- optimal in the sens of maximizing the sensitivity magnitude for each respective parameter 

% Outputs of the script:
% (1) Cardinality plot, indicating # exp in input library for which a
% parameter is identifiable according to our collinearity and sens. mag.
% threshold
% (2) Parameter Sensitivity Magnitude ranking plot
% (3) results structure: final parameters to be identified (unsorted and
% sorted from high->low sens. magnitude) and corresponding max sens. input
% for that parameter
clearvars
close all
clc

%% define threshold
threshold = 0.7;
Np = 22;
removed=[];

% Setup params string cell containing all parameters of interest
params = {'$D_s^{^{\_}}$','$D_s^+$','$R_s^-$','$R_s^+$',...
    '$\sigma^{^{\_}}$','$\sigma^{+}$','$D_e(\cdot)$','$\varepsilon_e^{^{\_}}$','$\varepsilon_e^{sep}$','$\varepsilon_e^+$',...
    '$\kappa(\cdot)$','$t_{c}^{0}$','$\frac{d \ln f_{c/a}}{d \ln c_e}(\cdot)$','$k^{^{\_}}$','$k^+$','$R_f^{^{\_}}$','$R_f^+$','$c_{e_0}$',...
    '$E.D_s^-$','$E.D_s^+$','$E.k_n$','$E.k_p$'};

params_reduced = params; % for tracking parameters remaining as some are eliminated

%% Load Inputs & Parameters

% load parameters (p struct)
run /Users/ztakeo/Documents/GitHub/OED_ParamID/Codebase/param/params_NCA.m

% load param boundaries (bounds struct)
run /Users/ztakeo/Documents/GitHub/OED_ParamID/Codebase/param/params_bounds.m

% Directory location for sensitivity .mat files; *****make sure they only have
% the .mat files for the inputs in them
senspath = '/Users/ztakeo/Documents/GitHub/OED_ParamID/SensResults/Tmax45/';
% senspath = '/Users/ztakeo/Documents/GitHub/OED_ParamID/SensResults/Tmax60/';

% Set output filename
output_filename = 'max_sens_experiments_Tmax45.mat';
% output_filename = 'max_sens_experiments_Tmax60.mat';

% Load sens files
sens_files = dir(senspath);
% on Mac, use line below to ignore '.','..', and '.DS_Store' files that
% are loaded into r 
sens_files=sens_files(~ismember({sens_files.name},{'.','..','.DS_Store'}));
num_exp = length(sens_files);

%% extract info from experiments

%%% Run the below function the first time after recalculating
%%% sensitivities. find_admissable_experiment_sets then saves the results
%%% into exp_info.mat, which can be loaded from thereon.

% [STSnorm,Sens_Mag,sens,NT_mat,exp_ind] = find_admissable_experiment_sets(p,params,removed,Np,senspath,num_exp,sens_files,bounds);

load('/Users/ztakeo/Documents/GitHub/OED_ParamID/SensResults/exp_info_Tmax45.mat','STSnorm','Sens_Mag','sens','NT_mat','exp_ind');
% load('/Users/ztakeo/Documents/GitHub/OED_ParamID/SensResults/exp_info_Tmax60.mat','STSnorm','Sens_Mag','sens','NT_mat','exp_ind');

if num_exp ~= length(sens)
    error('Check that senspath and exp_info.mat file are consistent i.e. have same # of experiments')
end
%% find parameters removed from clustering
removed = findclusterparams(STSnorm,Sens_Mag, threshold,Np,num_exp);

disp('Removed the following parameters during collinearity clustering:');
for pp = 1:length(removed)
    fprintf('%s \n', params{removed(pp)});
end
fprintf('\n');

% params_reduced(removed) = []; % for tracking purposes, remove clustered params from params cell



%% Remove parameters according to sensitivity threshold 

%update values will remove the params below the sensitivity threshold
[A_new,~] = updatevalues(sens, removed, NT_mat, Np, sens_files);

%%% Plot Cardinality and Orth Sensitivity after just removing collinearity
%%% clustered params (debugging)
% plotcardinality(A_new, params,removed,Np);
% plotcsensorth(Sens_orth_new, params, removed, -6, Np);

%parameters
remain_p = remaining(Np, removed);
card = Cardinality(A_new);
toremove = find(card == 0);
removed = [removed, remain_p(toremove)];
disp('Removed the following parameters according to noise sensitivity threshold');
for zz = 1:length(toremove)
    fprintf('%s \n', params{remain_p(toremove(zz))});
end

%% re-run gram-schmidt
[A_new,Sens_orth_new] = updatevalues(sens, removed, NT_mat,Np,sens_files);

%% determine paramters remaining
Np_final = Np - length(removed); 
params_final_idx = remaining(Np,removed);
params_remaining = params(params_final_idx)';

%% RePlot Cardinality & Orth Sensitivity for Remaining Params
plotcardinality(A_new, params,removed,Np);

param_sens_ranking = plotcsensorth(Sens_orth_new, params, removed, -5,Np);
param_sens_ranking = fliplr(param_sens_ranking)'; % note that the sens. ranking comes out lowest to highest; use fliplr to flip

%% Find Experiment that Maximizes Sensitivity for each param
% determine final number of parameters left; should be same as # rows in A_new and sens_orth_new
max_sens = zeros(size(Np_final));
max_exp_idx = zeros(size(Np_final));

% Iterate through orthogonalized sensitivity matrix rows (i.e. params)
for ii = 1:Np_final
    [max_sens(ii),max_exp_idx(ii)] = max(Sens_orth_new(ii,:));
end

% max_exp_idx gives the index (columwise from 1:1587) in the Sens_orth_new matrix corresponding to
% the input that maximizes sensitivity for that parameter; however, the
% column index in this matrix != the actual experiment number (i.e. exp. #10 isn't in the 10th spot of the experiment array. 
% To get this, need to grab the experiment numbers in exp_ind corresponding to this indice value (max_exp_idx). 
% This discrepancy is because some inputs are thrown out during the process
% due to constraint violations etc, making the experiment set
% non-contiguous. Additionally, filenames are not read in from lowest to
% highest or vice versa
max_exp_num = exp_ind(max_exp_idx)';

%% Sort parameters and their experiments according to ortho sens. ranking
params_sorted = params(params_final_idx(param_sens_ranking))';
max_exp_sorted = max_exp_num(param_sens_ranking);

%%% Check that each parameter in the group is above the sens threshold for
%%% each input
A_sort = A_new(param_sens_ranking,:);
group1_size = 7;
group2_size = 7;

fprintf('\n')
disp('Check that all parameters in Group 1 are identifiable for each input in that group.')
for kk = 1:group1_size
    for jj = 1:group1_size
        if (sum(A_sort(kk,:) == max_exp_sorted(jj)) == 0)
            fprintf('%s not sensitive to %i. Remove experiment. \n',params_sorted{kk}, max_exp_sorted(jj))
        end
    end
end
disp('Complete')
fprintf('\n')
disp('Check that all parameters in Group 2 are identifiable for each input in that group.')
for mm = 1:group2_size
    for pp = 1:group2_size
        if (sum(A_sort(mm,:) == max_exp_sorted(pp)) == 0)
            fprintf('%s not sensitive to %i. Remove experiment. \n',params_sorted{mm}, max_exp_sorted(jj))
        end
    end
end
disp('Complete')
%% Compile results
results.params_remaining = params_remaining;
results.max_exp_num = max_exp_num;
results.params_sorted = params_sorted;
results.max_exp_num_sorted = max_exp_sorted;
% save(output_filename,'results');
