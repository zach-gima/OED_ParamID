%% OED-CVX Parameter Estimation Function
% Adapted from Saheong Park Exp_param_ID_normalized_saehong_G1G2G3G4.m
% By Zach Gima, 2018-2-22

function [park0, paramID_out, Iter] = Param_ID(p,bounds,sel_k,selection_vector,theta_0,Inputs,filename_output,SCD_options)
%%
%p = parameter vector
%bounds = max-min bounds on parameters
%sel_k = a vector that contains the indicies corresponding to the vector of
%parameters
%selection vector: a vector of length(Num of paremeters total) that is 1 in
%the correspoinding index if that parameter is identified and zero if not

%% Setup Parameter Variables

% Setup vectors related to parameter sensitivity selection
Selected_params = theta_0(sel_k);  % Goes from 21x1 vector to 18x1 (in all params selected scenario)


%% Load Input and Set Variables

Current_exp = Inputs.Current_exp;
Time_exp = Inputs.Time_exp; %% room for memory improvement
Voltage_exp = Inputs.V_LM_CELL;
T_amb = Inputs.T_amb_sim; % note, comes in celcius
exp_num = Inputs.exp_num;

num_inputs = length(Voltage_exp);
v_dat = cell2mat(Voltage_exp); % Truth/experimental data
total_NT  = size(v_dat,1);

% In experimental ID, Rc needs to be identified for each experiment
% (Rc_tune) and the Rc value must be attached to each profile
if isfield(Inputs,'Rc')
    Rc = Inputs.Rc;
else %M2M case
    Rc = cell(length(Current_exp),1); % Create a cell num_exp x 1
    Rc(:,1) = {p.R_c}; % For M2M case, just use nominal Rc value
end


%% Normalization

normalized_theta_bar = origin_to_norm('param',Selected_params,bounds,selection_vector);

Min_theta_bar = zeros(25,1); % Year 2: Total # possible parameters = 25 (22 + 3 eq. struct)
Max_theta_bar = ones(25,1); % Year 2: Total # possible parameters = 25 (22 + 3 eq. struct)
normalize_params_min = Min_theta_bar(sel_k);
normalize_params_max = Max_theta_bar(sel_k);

normalized_sens_bar = origin_to_norm('sens',Selected_params,bounds,selection_vector);

%% SCM
groupsize = 1;
theta = theta_0;
Jac = 0;
exit_logic = false;

while exit_logic == false
    % Reset alpha 
    alpha = 1e-3;

    %Sample with replacement
    rand_idx = sel_k(randsample(sum(selection_vector),groupsize))
    e_idx = zeros(size(selection_vector));
    e_idx(rand_idx) = 1;
    
    %Sample without replacement
    %sel_k = ...randsample
    
    btrk = true;
    while btrk
        try
            V_CELL = cell(num_inputs,1);
            S_CELL = cell(num_inputs,1);
%             parfor idx = 1:num_inputs
            for idx = 1:num_inputs
                [V_CELL{idx}, ~, S_CELL{idx}] = DFN_sim_casadi(p,...
                    exp_num{idx},Current_exp{idx}, Time_exp{idx}, ...
                    Voltage_exp{idx}, T_amb{idx}, e_idx, theta, 1,Rc{idx});
            end
            v_sim = cell2mat(V_CELL);
            Sens = cell2mat(S_CELL);
            Jac = bsxfun(@times,normalized_sens_bar',Sens);
            
            % Update / increase alpha?            
            delta_theta = - alpha*(Jac')*(v_dat - v_sim);
            theta_prev = theta;
            theta = theta + delta_theta;
            
            %Check that the normalization is correct - Dylan/Zach
            btrk = false;
            
        catch % logic for when casadi fails
            % Propose a smaller parameter step
            alpha = alpha/100;
            delta_theta = - alpha*(Jac')*(v_dat - v_sim);
            theta = theta_prev + delta_theta;
            
        end
        
        
    end
    
    % Check Exit Conditions
    [exit_logic] = check_ec(v_dat,v_sim,delta_theta,Iter,SCD_options);
    
    % Save ParamID results every iteration
    %     paramID_out.Time_exp = Time_exp;
    %     paramID_out.Current_exp = Current_exp;
    %     paramID_out.Voltage_exp = Voltage_exp;
    %     paramID_out.save_chi_sq = chi_sq;
    %     paramID_out.save_L2norm = save_L2norm;
    %     paramID_out.save_RMSE = save_RMSE;
    %     paramID_out.save_delta_matrix = save_delta_matrix;
    %     paramID_out.save_param_nmz = save_param_nmz;
    %     paramID_out.save_param_org = save_param_org;
    %     paramID_out.save_y_minus_yfit = save_y_minus_yfit;
    %     paramID_out.save_y_sim = save_y_sim;
    %     paramID_out.save_param_exit = save_param_exit;
    %     paramID_out.save_chi_sq_exit = save_chi_sq_exit;
    %     paramID_out.save_chi_sq_AbsTol = save_chi_sq_AbsTol;
    %     paramID_out.save_lambda_matrix = save_lambda_matrix;
    %     paramID_out.y_dat = y_dat;
    %     paramID_out.LM_logic = LM_logic;
    
    %     save(strcat(LM_filename_output,num2str(LM_Iter),'.mat'),'park0','paramID_out','LM_Iter');
end

disp('After estimation')
fprintf('Final Parameter Values: %1.5e \n',park0)

%% Concatenate Outputs & Save Results
% paramID_out.exp_num = exp_num;
% paramID_out.Time_exp = Time_exp;
% paramID_out.Current_exp = Current_exp;
% paramID_out.Voltage_exp = Voltage_exp;
% paramID_out.save_chi_sq = chi_sq;
% paramID_out.save_L2norm = save_L2norm;
% paramID_out.save_RMSE = save_RMSE;
% paramID_out.save_delta_matrix = save_delta_matrix;
% paramID_out.save_param_nmz = save_param_nmz;
% paramID_out.save_param_org = save_param_org;
% paramID_out.save_y_minus_yfit = save_y_minus_yfit;
% paramID_out.save_y_sim = save_y_sim;
% paramID_out.save_param_exit = save_param_exit;
% paramID_out.save_chi_sq_exit = save_chi_sq_exit;
% paramID_out.save_chi_sq_AbsTol = save_chi_sq_AbsTol;
% paramID_out.save_lambda_matrix = save_lambda_matrix;
% paramID_out.y_dat = y_dat;
% paramID_out.LM_logic = LM_logic;

end