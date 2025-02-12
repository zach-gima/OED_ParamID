%% OED-CVX Parameter Estimation Main File
% Used to call ParamID function for OED-CVX param. estimation
% By Zach Gima, 2018-3-18
clc
clearvars
close all

datetime_initial = datetime('now','TimeZone','America/Los_Angeles');

%% User Input -- could also setup to just run both approaches and all I.C. variants;

% Load Nominal Parameter Set, Bounds, & True params
run param/params_NCA % loads p struct
run param/params_bounds % loads bounds struct (has fields bounds.min and bounds.max)
run param/params_truth % loads truth_param array

% Set Levenberg-Marquardt Conditions
param_exit_thresh = 1e-3; % 1e-6 default values used in matlab (called StepTol)
chi_sq_exit_thresh = 1e-6; % 1e-6 default values used in matlab (called FuncTol)
chi_sq_Abs_exit_thresh = 1e-4; % Absolute cost function tolerance (Func Absolute Tol)

LM_options.exit_cond = [param_exit_thresh, chi_sq_exit_thresh, chi_sq_Abs_exit_thresh];
LM_options.maxIter = 20;
LM_options.ctrl_lambda = 100; %1e-2; % initial lambda value (design variable); smaller = more optimistic and bigger initial steps -- SHP used 100 in yr 1

%%%%%%%%%%%%%%%   ParamID baseline (Uncomment to select)   %%%%%%%%%%%%%%%
% % Baseline A: Full Parameter Set (1 Group)
% baseline = {'full_'};
% num_groups = 1; % Number of parameter groups

% % Baseline B: Collinearity Only (1 Group)
% baseline = {'collinearity_'};
% num_groups = 1; % Number of parameter groups

% Baseline C: Collinearity + Sensitivity (2 Groups)
baseline = {'OED_'};
num_groups = 2; % Number of parameter groups

% % Experimental ParamID (Pre-Q Inclusion)
% baseline = {'OED_EXP_'};
% num_groups = 3; % Number of parameter groups G1(easy) G1&G2(easy) G1&G2(all)

%%%%%%%%%%%%%%%   Parameter Initial Conditions (Uncomment to select)   %%%%%%%%%%%%%%%
run param/params_nominal
init_cond = 'nom';
theta_0 = Nominal_param;

%%%%%%%%%%%%%%%  File I/O (Set once)   %%%%%%%%%%%%%%%
% Input subfolder
input_folder = strcat('InputLibrary/MaxSensInputs/OED/');
% input_folder = strcat('InputLibrary/MaxSensInputs/BaselineA/');
% input_folder = strcat('InputLibrary/MaxSensInputs/BaselineB/');
% input_folder = strcat('InputLibrary/ValidationCycles/');
% input_folder = strcat('InputLibrary/Experimental/');

% input_folder = strcat('InputLibrary/MaxSensInputs/plus50/');
% input_folder = strcat('InputLibrary/MaxSensInputs/minus50/');

% Output subfolder
date_txt = strrep(datestr(datetime_initial), ':', '_');
output_folder = strcat('/Users/ztakeo/Documents/GitHub/OED_ParamID/ID_results/',date_txt,'/');
% output_folder = strcat('C:/Users/Zach/Box Sync/HPC/HPC1/',date_txt,'/'); %HPC-1 Path
% output_folder = strcat('C:/Users/zgima/Box Sync/HPC/HPC2/',date_txt,'/'); %HPC-2 Path
% output_folder = strcat('/global/home/users/ztakeo/output/',date_txt,'/'); %Savio Path

mkdir(output_folder); %create new subfolder with current date in output_folder

%%% init_ParamID: initialize background stuff (variables, file i/o etc) based on the ParamID baseline and I.C.'s 
[filename_input_vector,filename_output_vector,selection_vector,ci_select,ci_input_vector] = init_ParamID(baseline,init_cond,num_groups,input_folder,output_folder);

%% Display Simulation Info

%try/catch structure used to send email alert if program exits w/ error
% For saving errors:
error_filename = strcat(output_folder,'sim_log.txt');
diary(error_filename)

datetime_initial
fprintf('Initial Conditions: %s \n',init_cond);
fprintf('Baseline: %s \n \n',baseline{1});
fprintf('Number of Groups: %i \n \n',num_groups);
disp('Levenberg-Marquardt Params')
fprintf('Initial Lambda: %5.2e \n',LM_options.ctrl_lambda);
fprintf('Param Convergence Exit Condition: %5.2e \n',LM_options.exit_cond(1));
fprintf('Cost Function Rel. Convergence Exit Condition: %5.2e \n',LM_options.exit_cond(2));
fprintf('Cost Function Abs. Convergence Exit Condition: %5.2e \n',LM_options.exit_cond(3));

fprintf('Max Iterations: %i \n \n',LM_options.maxIter);

%% Perturbation Analysis %%%%%%%%%
% %%% Comment out this section to not run
% %%% All true params except group X; perturb group X 10% away from
% %%% true values. 
% truth_param = Nominal_param;
% 
% % change this to indicate how many groups are being perturbed; e.g. for G2G1, num_groups = 2
% num_perturbedgroups = 2;  % for running V_sim_debug or perturbation analysis
% 
% % Perturb parameters of interest all at the beginning
% % perturb_index = find(selection_vector(:,1)); % G1
% perturb_index = find(selection_vector(:,end)); %G1 & G2
% 
% perturb_factor = 1.05;
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
%         theta_0(ii) = bounds.max(ii); % set param value to bound
% %         bounds.min(ii) = theta_0(ii); % set bound to param value
%         fprintf('Parameter %i perturbed past upper bound \n', ii);
%     end
%     
%     % Initial Value < Lower Bound
%     if theta_0(ii) < bounds.min(ii)
%         theta_0(ii) = bounds.min(ii); % set param value to bound
% %         bounds.min(ii) = theta_0(ii); % set bound to param value
%         fprintf('Parameter %i perturbed past lower bound \n', ii);
%     end
%     
% end

%% Call ParamID function

% Note: cannot parallelize this for loop, because cumulative baseline
% depends on G1 being ID'ed first, then G2 and so on
% Note: To use for Savio, probably need a separate main script with the
% Savio header

%initialize vectors for storing metrics and other data
datetime_paramID = cell(num_groups,1);
t_paramID = 0;
ci95_full = zeros(25,1); %vector for storing confidence interval for params
rmse_final = [];
iter_history = [];
theta_0_true = theta_0;% save the very 1st initial parameter guess for plotting purposes (param_table_plotter)

try
    for jj = 1:num_groups
        fprintf('Beginning Baseline %s, Group %s \n\n\n',baseline{1}(1:end-1),num2str(jj));
        
%         %% Change lambda depending on the group being identified (see note earlier about why)
%         LM_options.ctrl_lambda = ctrl_lambda(jj);
        
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
        t_paramID = vertcat(t_paramID,toc+t_paramID(end)); %measure program execution time
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
        'truth_param','theta_0_true','LM_options','sel_k','selection_vector','bounds','ci95_full','t_paramID','output_folder');
        
        % matlabmail(recipient,subject,message,attachments)
        matlabmail('ztakeo@berkeley.edu','Parameter ID complete','',[]);
        
     	%% Calculate Confidence Intervals at the last stage and Plot results
        disp('Calculate confidence intervals for identified parameter values \n')

        % Calculate local sensitivity and final RMSE
        tic
        [J_LM,ci95,sigma_y,covar_p,alg_states] ...
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
        'iter_history','output_folder','sel_k','selection_vector','truth_param',...
        'theta_0_true','LM_options','J_LM','bounds','alg_states');

        matlabmail('ztakeo@berkeley.edu','Conf. Interval complete','',[]);
    end
    
    %% Plot Results: plot figure of truth and estimated values w/ C.I.'s
%     Param_ID_plot(truth_param,theta_0_true,sel_k,paramID_out,ci95_full,t_paramID,rmse_final,output_folder,LM_options,bounds,alg_states,selection_vector)    

catch e %e is an MException struct
    % An error will put you here.
    errorMessage = sprintf('%s',getReport( e, 'extended', 'hyperlinks', 'on' ))
    
%     % Save ParamID results even when erorr occurs
%     save(filename_output,'park0','paramID_out','LM_Iter','ci95_full');
    
    diary off %stop logging command window
    
    matlabmail('ztakeo@berkeley.edu','Error encountered',errorMessage,error_filename);
end

