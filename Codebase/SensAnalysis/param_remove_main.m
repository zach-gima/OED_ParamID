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

%% Set Group Sizes
group1_size = 6;
group2_size = 7;

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
% Set output folder
output_folder = '/Users/ztakeo/Documents/GitHub/OED_ParamID/Codebase/Plots/SensAnalysis/Tmax60/';
% output_folder = '/Users/ztakeo/Documents/GitHub/OED_ParamID/Codebase/Plots/Perturb/SensAnalysis/minus50/';
% output_folder = '/Users/ztakeo/Documents/GitHub/OED_ParamID/Codebase/Plots/Perturb/SensAnalysis/plus50/';

% load parameters (p struct)
run /Users/ztakeo/Documents/GitHub/OED_ParamID/Codebase/param/params_NCA.m

% load param boundaries (bounds struct)
run /Users/ztakeo/Documents/GitHub/OED_ParamID/Codebase/param/params_bounds.m

% Directory location for sensitivity .mat files; *****make sure they only have
% the .mat files for the inputs in them
% senspath = '/Users/ztakeo/Documents/GitHub/OED_ParamID/SensResults/Tmax45/';
senspath = '/Users/ztakeo/Documents/GitHub/OED_ParamID/SensResults/Tmax60/';

% Set output filename
% output_filename = 'max_sens_experiments_Tmax45.mat';
output_filename = 'max_sens_experiments_Tmax60.mat';
% output_filename = 'max_sens_experiments_minus50.mat';
% output_filename = 'max_sens_experiments_plus50.mat';

% Load sens files
sens_files = dir(senspath);
% on Mac, use line below to ignore '.','..', and '.DS_Store' files that
% are loaded into r 
sens_files=sens_files(~ismember({sens_files.name},{'.','..','.DS_Store'}));
num_exp = length(sens_files);


%% Perturbation Analysis %%%%%%%%%
% %%% Comment out this section to not run
% %%% All true params except group X; perturb group X 10% away from
% %%% true values. 
% run /Users/ztakeo/Documents/GitHub/OED_ParamID/Codebase/param/params_nominal.m
% theta_0 = Nominal_param;
% 
% % Load feasible inputs for perturbed parameter sets
% % perturb_path = '/Users/ztakeo/Documents/GitHub/OED_ParamID/SensResults_Perturb/minus50.mat';
% perturb_path = '/Users/ztakeo/Documents/GitHub/OED_ParamID/SensResults_Perturb/plus50.mat';
% load(perturb_path);
% perturb_admissable = arr;
% 
% % change this to indicate how many groups are being perturbed; e.g. for G2G1, num_groups = 2
% num_perturbedgroups = 2;  % for running V_sim_debug or perturbation analysis
% 
% % Perturb parameters of interest all at the beginning
% % perturb_index = find(selection_vector(:,1)); % G1
% selection_vector(:,2) = [1;1;1;1;0;0;0;0;1;0;1;1;0;0;1;0;0;1;0;0;1;1;1;1;0]; %G2
% perturb_index = find(selection_vector(:,2)); %G1 & G2
% 
% % perturb_factor = 0.5; % minus50
% perturb_factor = 1.5; % plus50
% theta_0(perturb_index) = perturb_factor*theta_0(perturb_index);
% 
% % Display Analysis Info
% fprintf('Perturbation Analysis Experiment \n');
% fprintf('Groups 1->%i \n',num_perturbedgroups); 
% fprintf('Number of parameters: %i \n',length(perturb_index));
% fprintf('Perturb Factor = %5.2f \n \n', perturb_factor);
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
% 
% % Check whether the perturbed value exceeds a pre-set bound (params_bounds)
% % If so, then replace the bound with the initial value
% for ii = 1:length(theta_0)  
%     
%     % Initial Value > Upper Bound
%     if theta_0(ii) > bounds.max(ii)
%         bounds.max(ii) = theta_0(ii);
%         fprintf('Parameter %i perturbed past upper bound \n', ii);
%     end
%     
%     % Initial Value  Lower Bound
%     if theta_0(ii) < bounds.min(ii)
%         bounds.min(ii) = theta_0(ii);
%         fprintf('Parameter %i perturbed past lower bound \n', ii);
%     end
%     
% end

%% extract info from experiments

%%% Run the below function the first time after recalculating
%%% sensitivities. find_admissable_experiment_sets then saves the results
%%% into exp_info.mat, which can be loaded from thereon.

% [STSnorm,Sens_Mag,sens,NT_mat,exp_ind] = find_admissable_experiment_sets(p,params,removed,Np,senspath,num_exp,sens_files,bounds,output_folder);

% load('/Users/ztakeo/Documents/GitHub/OED_ParamID/SensResults/exp_info_Tmax45.mat','STSnorm','Sens_Mag','sens','NT_mat','exp_ind');
load('/Users/ztakeo/Documents/GitHub/OED_ParamID/SensResults/exp_info_Tmax60.mat','STSnorm','Sens_Mag','sens','NT_mat','exp_ind');

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
fig = plotcardinality(A_new, params,removed,Np);
savefig(fig,strcat(output_folder,'cardinality_final_Tmax60'));
saveas(fig,strcat(output_folder,'cardinality_final_Tmax60.png'));
% savefig(fig,strcat(output_folder,'cardinality_final_minus50'));
% saveas(fig,strcat(output_folder,'cardinality_final_minus50.png'));
% savefig(fig,strcat(output_folder,'cardinality_final_plus50'));
% saveas(fig,strcat(output_folder,'cardinality_final_plus50.png'));

[fig,param_sens_ranking] = plotcsensorth(Sens_orth_new, params, removed, -5,Np);
param_sens_ranking = fliplr(param_sens_ranking)'; % note that the sens. ranking comes out lowest to highest; use fliplr to flip
savefig(fig,strcat(output_folder,'orthsens_final_Tmax60'));
saveas(fig,strcat(output_folder,'orthsens_final_Tmax60.png'))
% savefig(fig,strcat(output_folder,'orthsens_final_minus50'));
% saveas(fig,strcat(output_folder,'orthsens_final_minus50.png'))
% savefig(fig,strcat(output_folder,'orthsens_final_plus50'));
% saveas(fig,strcat(output_folder,'orthsens_final_plus50.png'))

%% Find Experiment that Maximizes Sensitivity for each param
% determine final number of parameters left; should be same as # rows in A_new and sens_orth_new

%%%% Iterate through orthogonalized sensitivity matrix rows (i.e. params)
% no_sensitivity{1} = [];% inputs for which some or multiple group 1 params are not sensitive
% no_sensitivity{2} = [930 875 245 21 580 221 83 650 858 860 251 270 10 239 937 861 276 287 657 1077 617 215 664 593 867 32 ...
%     294 47 944 596 601];% inputs for which some or multiple group 2 params are not sensitive 

[Sens_orth_sort,sort_idx] = sort(Sens_orth_new,2,'descend');
A_sort = A_new(param_sens_ranking,:); % sort the A matrix according to the param sens. ranking

max_exp_idx = zeros(Np_final,1);
max_sens = zeros(Np_final,1);

for ii = 1:Np_final
    %%%% Routine for choosing max sensitivity experiment, not considering
    %%%% whether rest of the parameters in the group are sensitive to that
    %%%% input
    [max_sens(ii),max_exp_idx(ii)] = max(Sens_orth_new(ii,:));
    
    %%%% Routine considering whether rest of parameters in a group are
    %%%% sensitive to each input
%     for jj = 1:length(Sens_orth_sort) 
%         % Check Group 1
%         if ii <=group1_size
%             sens_check = 0;
%             for zz = 1:group1_size
%                 sens_check = sens_check + sum(A_sort(zz,:) == exp_ind(sort_idx(ii,jj)));
%             end
%             if sens_check == group1_size %~ismember(exp_ind(sort_idx(ii,jj)),no_sensitivity{1})
%                 max_sens(ii) = Sens_orth_sort(ii,jj);
%                 max_exp_idx(ii) = sort_idx(ii,jj);
%                 break;
%             end
%         % Check Group 1
%         else
%             sens_check = 0;
%             for zz = (group1_size+1):Np_final
%                 sens_check = sens_check + sum(A_sort(zz,:) == exp_ind(sort_idx(ii,jj)));
%             end
%             if sens_check == group2_size %~ismember(exp_ind(sort_idx(ii,jj)),no_sensitivity{2})
%                 max_sens(ii) = Sens_orth_sort(ii,jj);
%                 max_exp_idx(ii) = sort_idx(ii,jj);
%                 break;
%             end
%         end
%     end

end

% %%%%%% Method when using perturbed parameter set 
%%% Checking that the selected input is a member of the inputs for which
%%% the perturbed parameter set can actually run; the variable
%%% perturb_admissable contains the experiment numbers that are feasible
% [Sens_orth_sort,sort_idx] = sort(Sens_orth_new,2,'descend');
% 
% for ii = 1:Np_final
%     for jj = 1:length(Sens_orth_sort)
%         if ismember(exp_ind(sort_idx(ii,jj)),perturb_admissable)
%            max_sens(ii) = Sens_orth_sort(ii,jj);
%            max_exp_idx(ii) = sort_idx(ii,jj);
%            break;
%         end
%     end
% end

% max_exp_idx gives the index (columwise from 1:num_exp) in the Sens_orth_new matrix corresponding to
% the input that maximizes sensitivity for that parameter; however, the
% column index in this matrix != the actual experiment number (i.e. exp. #10 isn't in the 10th spot of the experiment array. 
% To get this, need to grab the experiment numbers in exp_ind corresponding to this indice value (max_exp_idx). 
% This discrepancy is because some inputs are thrown out during the process
% due to constraint violations etc, making the experiment set
% non-contiguous. Additionally, filenames are not read in from lowest to
% highest or vice versa
max_exp_sorted = exp_ind(max_exp_idx)';

%% Sort parameters and their experiments according to ortho sens. ranking
params_sorted = params(params_final_idx(param_sens_ranking))';
% max_exp_sorted = max_exp_num(param_sens_ranking);

%check that group sizes set add up to total number of params to be
%identified
if group1_size + group2_size ~= length(params_sorted)
    error('Check group sizes and make consistent with number of parameters remaining after sensitivity and clustering analysis');
end

%%% Check that each parameter in the group is above the sens threshold for
%%% each input.
%%% Routine: A_sort is a [num_param x num_exp] sized matrix, where the rows are ordered according to sensitivity (highest sens. param is 1st row, 
%%% 2nd highest is 2nd and so on. Each row and column pair indicates whether the parameter is sensitive for a given
%%% experiment (index listed in the array value). To check that each
%%% parameter is sensitive to each experiment in its group, we loop through
%%% each row of A_sort and take the sum of the logical expression checking
%%% whether the chosen experiment index shows up. For example, if Exp. 488
%%% is chosen, for each parameter in Group 1 to be sensitive to that input,
%%% the expression sum(A_sort(kk,:) == max_exp_sorted(jj)) should equal 1.
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

% For group 2, be careful with setting the for loop indices. We are now
% looping through the rows of A_sort corresponding to the group 2 params,
% which start at row group1_size + 1 and end at the number of params
% remaining. e.g. if G1 = 6, G2 = 7, then we start at 7 and end at 13.
disp('Check that all parameters in Group 2 are identifiable for each input in that group.')
for mm = (group1_size+1):Np_final
    for pp = (group1_size+1):Np_final
        if (sum(A_sort(mm,:) == max_exp_sorted(pp)) == 0)
            fprintf('%s not sensitive to %i. Remove experiment. \n',params_sorted{mm}, max_exp_sorted(pp))
        end
    end
end
disp('Complete')

%% Compile results
results.params_remaining = params_remaining;
% results.max_exp_num = max_exp_num;
results.params_sorted = params_sorted;
results.max_exp_num_sorted = max_exp_sorted;
save(output_filename,'results');
