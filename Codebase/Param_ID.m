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

%% Load Input and Set Variables

Current_exp = Inputs.Current_exp;
Time_exp = Inputs.Time_exp; %% room for memory improvement
Voltage_exp = Inputs.V_LM_CELL;
T_amb = Inputs.T_amb_sim; % note, comes in celcius
exp_num = Inputs.exp_num;

% num_inputs = length(Voltage_exp);
num_inputs = 1;
% v_dat = cell2mat(Voltage_exp); % Truth/experimental data
v_dat = Voltage_exp{2}(1:end);
total_NT  = size(v_dat,1);

% In experimental ID, Rc needs to be identified for each experiment
% (Rc_tune) and the Rc value must be attached to each profile
if isfield(Inputs,'Rc')
    Rc = Inputs.Rc;
else %M2M case
    Rc = cell(length(Current_exp),1); % Create a cell num_exp x 1
    Rc(:,1) = {p.R_c}; % For M2M case, just use nominal Rc value
end

%% SCM
groupsize = 1;
theta = theta_0;
Jac = 0;
exit_logic = false;

% DEBUG
plot(v_dat,'LineWidth',2,'Color','k');
hold on

np = sum(selection_vector);
delta_theta_history = ones(np,1)*100;

while exit_logic == false
    % Reset alpha
    alpha = 1e5;
    
    %Sample with replacement
    rand_idx25 = sel_k(randsample(np,groupsize))
    %     rand_idx25 = 25;
    rand_idx22 = twenty2to25(rand_idx25,'522');
    e_idx = zeros(size(selection_vector));
    e_idx(rand_idx25) = 1;
    
    %Sample without replacement
    %sel_k = ...randsample
    
    % Pick random coordinate, optimize cost function via approx. line search along that
    % coordinate
    btrk = true;
    
    p = update_p(p,theta);
    
    V_CELL = cell(num_inputs,1);
    S_CELL = cell(num_inputs,1);
    
    %             parfor idx = 1:num_inputs
    for idx = 1:num_inputs
        [V_CELL{idx}, ~, S_CELL{idx}] = DFN_sim_casadi(p,...
            exp_num{idx+1},Current_exp{idx+1}(1:end), Time_exp{idx+1}(1:end), ...
            Voltage_exp{idx+1}(1:end), T_amb{idx+1}, e_idx, theta, 1,Rc{idx+1});
    end
    
    v_sim = cell2mat(V_CELL);
    Sens = cell2mat(S_CELL);
    costprev = norm(v_dat - v_sim,2);
    
    % Normalize Sensitivity
    Selected_params = theta(sel_k);  % Goes from 25x1 vector to 22x1 (in all params selected scenario)
    normalized_sens_bar = origin_to_norm('sens',Selected_params,bounds,selection_vector);
    Jac = bsxfun(@times,normalized_sens_bar(rand_idx25),Sens);
    
    alpha =min(1e3 ,1/(Jac')*(v_dat - v_sim));
    
    theta_prev = theta; % NOTE: UN-NORMALIZED this is the theta from the previous successful iteration
    
    %% Line search
    while btrk
        try
            while true
                %             update parameter, just simulate to see if cost improves,
                %             if not reduce alpha,
                %             if so, break
                
                % DEBUG
                plot(v_sim)
                drawnow
                
%                 % Normalize Sensitivity
%                 Selected_params = theta(sel_k);  % Goes from 25x1 vector to 22x1 (in all params selected scenario)
%                 normalized_sens_bar = origin_to_norm('sens',Selected_params,bounds,selection_vector);
%                 Jac = bsxfun(@times,normalized_sens_bar(rand_idx25),Sens);
                
                % Update / increase alpha?
                delta_theta = alpha*(Jac')*(v_dat - v_sim); % NOTE: THIS UPDATE IS NORMALIZED
                delta_theta_history(rand_idx22) = delta_theta;
                
                % Parameter Normalization -- NOTE: NEEDS TESTING
                %theta_prev = theta; % NOTE: UN-NORMALIZED this is the theta from the previous successful iteration
                theta_norm = origin_to_norm('param',Selected_params,bounds,selection_vector);
                theta_norm(rand_idx22) = theta_norm(rand_idx22) + delta_theta;
                theta_norm(rand_idx22) = min(max(0,theta_norm(rand_idx22)),1); % min max routine prevents parameter value from violating bounds
                theta = norm_to_origin(theta_norm,bounds,selection_vector);
                
                %check if we should got to another parameter or update this one
                %again
                V_CELL = cell(num_inputs,1);
                S_CELL = cell(num_inputs,1);
                
                for idx = 2%1:num_inputs
                    [V_CELL{idx-1}] = DFN_sim_casadi(p,...
                        exp_num{idx},Current_exp{idx}(1:end), Time_exp{idx}(1:end), ...
                        Voltage_exp{idx}(1:end), T_amb{idx}, e_idx, theta, 0,Rc{idx});
                end
                v_new = cell2mat(V_CELL);
                costnew = norm(v_dat - v_new,2)
                
                %check if the cost got worse
                if costnew > costprev
                    %reset theta
                    theta = theta_prev;
                    %reduce step size
                    grad = (Jac')*(v_dat - v_sim);
                    alpha = min(alpha,1/abs(grad));
                    alpha = alpha/2;
                else
                    %It will reach this when a step improves the cost
                    btrk = false;
                    break
                end
            end
            
            
            
        catch e % logic for when casadi fails %TEST THIS
            % An error will put you here.
            errorMessage = sprintf('%s',getReport( e, 'extended', 'hyperlinks', 'on' ))
            
            % Propose a smaller parameter step
            alpha = alpha/2; % NOTE: NEEDS TESTING
            %delta_theta = alpha*(Jac')*(v_dat - v_sim); % NOTE: THIS UPDATE IS NORMALIZED
            %grad = (Jac')*(v_dat - v_sim);
            %alpha = min(alpha,1/abs(grad));
            
            % Parameter Normalization -- NOTE: NEEDS TESTING
            %theta_norm = origin_to_norm('param',Selected_params,bounds,selection_vector);
            %theta_norm(rand_idx22) = theta_prev(rand_idx22) + delta_theta;
            %theta_norm(rand_idx22) = min(max(0,theta_norm(rand_idx22)),1); % min max routine prevents parameter value from violating bounds
            %theta = norm_to_origin(theta_norm,bounds,selection_vector);
            theta = theta_prev;
        end
        
        
    end
    
    % Check Exit Conditions
    %     [exit_logic] = check_ec(v_dat,v_sim,delta_theta_history,Iter,SCD_options);
    
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
    cost = norm(v_dat - v_new,2);
    fprintf('Cost: %1.6f \n',cost);
    
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