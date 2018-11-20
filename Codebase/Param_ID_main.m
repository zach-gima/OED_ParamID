%% OED-CVX Parameter Estimation Main File
% Used to call ParamID function for OED-CVX param. estimation
% By Zach Gima, 2018-3-18
clc
clearvars
close all

datetime_initial = datetime('now','TimeZone','America/Los_Angeles')

%Load CasADi
% addpath('C:/Users/Zach/Documents/MATLAB/casadi_windows')
% import casadi.*

%% User Input -- could also setup to just run both approaches and all I.C. variants;

% Load Nominal Parameter Set, Bounds, & True params
run param/params_NCA % loads p struct
run param/params_bounds % loads bounds struct (has fields bounds.min and bounds.max)
run param/params_truth % loads truth_param array

% Set Levenberg-Marquardt Conditions
param_exit_thresh = 1e-6; % 1e-6 default values used in matlab (called StepTol)
chi_sq_exit_thresh = 1e-6; % 1e-6 default values used in matlab (called FuncTol)

LM_options.exit_cond = [param_exit_thresh, chi_sq_exit_thresh];
LM_options.maxIter = 40;

LM_options.ctrl_lambda = 1e-2; %100; % initial lambda value (design variable) -- SHP used 100 in yr 1

%%% Number of parameter groups
num_groups = 2;

%%%%%%%%%%%%%%%   ParamID Approach (Uncomment to select)   %%%%%%%%%%%%%%%

% (1) Cumulative
approach = {'G1_'; 'G2G1_'};
init_iter = 1; % start w/ G1 params

% (2) All-at-once (Use for testing with true initial parameters)
% approach = {'';'all_'};
% init_iter = num_groups; % identify all parameters together, so just 1 iteration needed 


%%%%%%%%%%%%%%%   Parameter Initial Conditions (Uncomment to select)   %%%%%%%%%%%%%%%
% (1) Theta_0 = True parameter values ID'ed by SHP in Bosch Year 1 (NCA)
%%% Note: will need to use sensitivity vector = all parameters for sens_lm
%%% to simulate voltage with all of the true params.

% init_cond = 'true';
% theta_0 = truth_param;

% (2) Theta_0 = Nominal_param from literature; Used in Bosch Year 1 (NCA)

run param/params_nominal
init_cond = 'nom';
theta_0 = Nominal_param;

% (3) Theta_0 = param lower bound
% Initialize close to but not at the bound

% init_cond = 'lb';
% theta_0 = param_min*1.25;

% (4) Theta_0 = param upper bound
% Initialize close to but not at the bound

% init_cond = 'ub';
% theta_0 = param_max*0.50;


%%%%%%%%%%%%%%%  File I/O (Set once)   %%%%%%%%%%%%%%%
% Input subfolder
input_folder = strcat('InputLibrary/MaxSensInputs/Tmax60/');

% Output subfolder
date_txt = strrep(datestr(datetime_initial), ':', '_');
output_folder = strcat('/Users/ztakeo/Documents/GitHub/OED_ParamID/ID_results/',date_txt,'/');

mkdir(output_folder); %create new subfolder with current date in output_folder
% output_folder = strcat(io_folder,'ID_results/',strrep(datestr(datetime_initial), ':', '_'),'/'); %rename output folder with newly created subfolder

%%% init_ParamID: initialize background stuff (variables, file i/o etc) based on the ParamID approach and I.C.'s 
[filename_input_vector,filename_output_vector,selection_vector,ci_select,ci_input_vector] = init_ParamID(approach,init_cond,input_folder,output_folder);

%% Display Simulation Info
fprintf('Initial Conditions: %s \n \n',init_cond);

disp('Levenberg-Marquardt Params')
fprintf('Initial Lambda: %5.2e \n',LM_options.ctrl_lambda);
fprintf('Param Convergence Exit Condition: %5.2e \n',LM_options.exit_cond(1));
fprintf('Cost Function Convergence Exit Condition: %5.2e \n',LM_options.exit_cond(2));
fprintf('Max Iterations: %i \n \n',LM_options.maxIter);

%% Perturbation Analysis %%%%%%%%%
%%% Comment out this section to not run
%%% All true params except group X; perturb group X 10% away from
%%% true values. 

% % change this to indicate how many groups are being perturbed; e.g. for G2G1, num_groups = 2
% num_groups = 3;  % for running V_sim_debug or perturbation analysis
% 
% % Perturb parameters of interest all at the beginning
% % perturb_index = [3;4]; % G1 only
% % perturb_index = [1;2;3;4;9;10;13;15]; % G1 & G2
% perturb_index = [1;2;3;4;9;10;12;13;15;16;18;19;21]; % G1, G2, G3
% 
% perturb_factor = 1.3;
% theta_0(perturb_index) = perturb_factor*theta_0(perturb_index);
% 
% % Display Analysis Info
% fprintf('Perturbation Analysis Experiment \n');
% fprintf('Groups 1->%i \n',num_groups); 
% fprintf('Number of parameters: %i \n',length(perturb_index));
% fprintf('Perturb Factor = %5.2f \n \n', perturb_factor);
% 
% % Overwrite the nominal values in the param struct w/ true values
% p.D_s_n0 = theta_0(1); %[G2]
% p.D_s_p0 = theta_0(2); %[G2]
% p.R_s_n = theta_0(3); %[G1]
% p.R_s_p = theta_0(4); %[G1]
% p.epsilon_s_n = theta_0(5); %  Equil. Struct
% p.epsilon_s_p = theta_0(6); %  Equil. Struct
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
% p.n_Li_s = theta_0(20); %  Equil. Struct
% p.c_e0 = theta_0(21);
% p.E.Dsn = theta_0(22);
% p.E.Dsp = theta_0(23);
% p.E.kn = theta_0(24);
% p.E.kp = theta_0(25);
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
%     % Initial Value < Lower Bound
%     if theta_0(ii) < bounds.min(ii)
%         bounds.min(ii) = theta_0(ii);
%         fprintf('Parameter %i perturbed past lower bound \n', ii);
%     end
%     
% end

%% Call ParamID function

% Note: cannot parallelize this for loop, because cumulative approach
% depends on G1 being ID'ed first, then G2 and so on
% Note: To use for Savio, probably need a separate main script with the
% Savio header

%try/catch structure used to send email alert if program exits w/ error
% For saving errors:
error_filename = strcat(output_folder,'error.txt');
diary(error_filename)

%initialize vectors for storing metrics and other data
datetime_paramID = cell(num_groups,1);
t_paramID = 0;
% ci95_full = zeros(21,1); %vector for storing confidence interval for params
ci95_full = zeros(25,1); %vector for storing confidence interval for params
rmse_final = [];
iter_history = [];
theta_0_true = theta_0;% save the very 1st initial parameter guess for plotting purposes (param_table_plotter)

try
    for jj = init_iter:num_groups
        fprintf('Beginning Stage %s \n\n\n',approach{jj});
        %%%%%%%% Debug %%%%%%%%%%%
%         filename_input_vector{1} = strcat(input_folder,'V_sim_debug_capiaglia'); %debugging input (much shorter)
%         selection_vector(:,1) = selection_vector(:,end); %selection vector = all params; use for debugging with true initial parms
        %%%%%%%%% Debug %%%%%%%%%%%%%%
        
        %% Load Group Specific Inputs
        filename_input = filename_input_vector{jj};
        filename_output = filename_output_vector{jj};
        Inputs = load(filename_input); %Current, Voltage, Time, T_amb     
        ci_inputs = load(ci_input_vector{jj}); %Inputs for each separate group -- used for C.I. calc.
        
        SensSelec = selection_vector(:,jj);
        sel_k = find(selection_vector(:,jj));

        %% Call function       
        tic %measure program execution time

        % Prev. had outputs [LM_Iter,paramID_out], which are now just
        % global variables. See note at beginning of script
        [park0, paramID_out, LM_Iter] = Param_ID(p,bounds,sel_k,selection_vector(:,jj),theta_0,Inputs,filename_output,LM_options); 

        % Save time it took to identify each parameter group
        datetime_paramID{jj} = datetime('now','TimeZone','local')
        t_paramID = vertcat(t_paramID,toc+t_paramID(jj-1)); %measure program execution time
        iter_history = vertcat(iter_history,LM_Iter);
       
        % Replace parameter values with newly ID'ed values
        theta_0(sel_k) = park0;
           
        % Save very initial RMSE
        if isempty(rmse_final)
            rmse_final = vertcat(rmse_final,paramID_out.save_RMSE(1));
        end
        
        % Calculate final RMSE
        disp('Calculate final RMSE for identified parameter values')
        
        rmse_final = vertcat(rmse_final,paramID_out.save_RMSE(end));
        fprintf('##### Final RMSE : %f ##### \n', rmse_final(end));
            
        % Save data & send email
        save(filename_output,'paramID_out','rmse_final','datetime_initial','datetime_paramID',...
        'truth_param','theta_0_true','LM_options','sel_k','bounds','ci95_full','t_paramID','output_folder');

        % consider deleting data from previous L-M iterations
        % LM_filename_output = strcat(filename_output(1:end-4),'_iter_');
        % delete LM_filename_output*
        
        % matlabmail(recipient,subject,message,attachments)
        matlabmail('ztakeo@berkeley.edu','Parameter ID complete','',[]);
        
     	%% Calculate Confidence Intervals at the last stage and Plot results
        %%%% Debug
        % Load park0 for group 4 if you had to exit manually
%         load('C:\Users\Zach\Documents\GitHub\DFN_CasADi\3_Journal_JPS\m2m_comp\ID_results\April 7\park0_G4.mat')
        %%%% Debug
        disp('Calculate confidence intervals for identified parameter values \n')

        % Calculate local sensitivity and final RMSE
        tic
        [J_LM,ci95,sigma_y,covar_p] ...
            = conf_interval(p, ci_inputs, SensSelec, park0);

        t_ci = toc; %store time for computing confidence interval

        %%% Store Confidence Intervals
        % After the 1st group, need to only replace the confidence intervals
        % calculated for the current group. This is necessary because
        % conf_interval calculates for all the identified parameters passed in.
        % For example, in G2G1 group, G1 and G2 conf intervals are calcualted,
        % but we've already calculated the correct G1 intervals the previous
        % iteration and don't want to overwrite them.
        ci95 = horzcat(ci95,sel_k); % append indices for each parameter; will be used below
        if jj > 1
            %iterate through the indices in ci_select
            for kk = 1:length(ci_select{jj})
                % iterate through the indices in sel_k, the 2nd column of ci95
                for mm = 1:size(ci95,1)
                    % when the sel_k index matches the ci_select index, then
                    % write that conf. interval value to ci95_full, which
                    % stores everything
                    if ci95(mm,2) == ci_select{jj}(kk)
                        ci95_full(ci_select{jj}(kk)) = ci95(mm,1); 
                        break;
                    end
                end
            end
        else
            ci95_full(ci_select{1}) = ci95(:,1);
        end

        datetime_final = datetime('now','TimeZone','local') % total program end time

        %fprintf('##### Final RMSE : %f ##### \n', rmse_final(jj));
        fprintf('The parameter estimate is %1.5e \n',park0)
        fprintf('plus/minus 95 %% confidence interval is %1.5e \n',ci95_full(sel_k))

        % Export Results
        save(filename_output,'paramID_out','rmse_final','ci95_full','sigma_y',...
        'covar_p','datetime_initial','datetime_final','t_paramID','t_ci',...
        'iter_history','output_folder','sel_k','truth_param','theta_0_true','LM_options','J_LM','bounds','alg_states');

        matlabmail('ztakeo@berkeley.edu','Parameter ID complete','',[]);

        %% Plot Results: plot figure of truth and estimated values w/ C.I.'s
%         Param_ID_plot(truth_param,theta_0_true,sel_k,paramID_out,ci95_full,t_paramID,rmse_final,output_folder,LM_options,bounds,alg_states)
        
    end

catch e %e is an MException struct
    % An error will put you here.
    errorMessage = sprintf('%s',getReport( e, 'extended', 'hyperlinks', 'on' ))
    
%     % Save ParamID results even when erorr occurs
%     save(filename_output,'park0','paramID_out','LM_Iter','ci95_full');
    
    diary off %stop logging command window
    
    matlabmail('ztakeo@berkeley.edu','Error encountered',errorMessage,error_filename);
end

