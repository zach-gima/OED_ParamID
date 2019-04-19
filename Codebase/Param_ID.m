%% OED-CVX Parameter Estimation Function
% Adapted from Saheong Park Exp_param_ID_normalized_saehong_G1G2G3G4.m
% By Zach Gima, 2018-2-22

function [park0, paramID_out, LM_Iter] = Param_ID(p,bounds,sel_k,selection_vector,theta_0,Inputs,filename_output,LM_options)
%%
%p = parameter vector
%bounds = max-min bounds on parameters
%sel_k = a vector that contains the indicies corresponding to the vector of
%parameters
%selection vector: a vector of length(Num of paremeters total) that is 1 in
%the correspoinding index if that parameter is identified and zero if not

%% Setup Parameter Variables
% Parse L-M Conditions
param_exit_thresh = LM_options.exit_cond(1);
chi_sq_exit_thresh = LM_options.exit_cond(2);
chi_sq_Abs_exit_thresh = LM_options.exit_cond(3);

maxIter = LM_options.maxIter;
ctrl_lambda = LM_options.ctrl_lambda;

% Setup vectors related to parameter sensitivity selection
SensSelec = selection_vector;
Selected_params = theta_0(sel_k);  % Goes from 21x1 vector to 18x1 (in all params selected scenario)
num_param = size(Selected_params,1); % number of param


%%%LM_logic = 0; % variable to track what actions L-M is taking (-1: decrease lambda; 0: recalc sensitivity; 1: increase lambda)

% Check things are imported properly
%%%if ( size(nonzeros(SensSelec),1) ~= size(Selected_params,1)); error('size is wrong'); end

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

%% Levenberg_Marquardt
norm_park0 = normalized_theta_bar;
groupsize = 1;
alpha = 1e-3;
theta = theta_0;


for LM_Iter=1:maxIter
    %Sample with replacement
    rand_idx = sel_k(randsample(sum(selection_vector),groupsize))
    e_idx = zero(size(selection_vector));
    e_idx(rand_idx) = 1;
    
    
    %Sample without replacement
    %sel_k = ...randsample
    
    btrk = true;
    while btrk
        try
            V_CELL = cell(num_inputs,1);
            S_CELL = cell(num_inputs,1);
            parfor idx = 1:num_inputs
                [V_CELL{idx}, ~, S_CELL{idx}] = DFN_sim_casadi(p,...
                    exp_num{idx},Current_exp{idx}, Time_exp{idx}, ...
                    Voltage_exp{idx}, T_amb{idx}, e_idx, theta, 1,Rc{idx});
            end
            v_sim = cell2mat(V_CELL);
            Sens = cell2mat(S_CELL);
            J_LM = bsxfun(@times,normalized_sens_bar',Sens);
            delta_theta = - alpha*(J_LM')*(v_dat - v_sim);
            theta_prev = theta;
            theta = theta + delta_theta;
            
            %Check that the normalization is correct
            
            btrk = false;
        catch
            %implemet logic for when casadi fails
        end
        
        
    end

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1. Param updates & Simulation Run
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    norm_park0 = norm_park0 + delta_theta;
    norm_park0 = min(max(normalize_params_min,norm_park0),normalize_params_max); % min max routine prevents parameter value from violating bounds
    save_param_nmz(:,LM_Iter)=norm_park0;
    save_param_org(:,LM_Iter)=norm_to_origin(norm_park0,bounds,selection_vector);
    
    fprintf('---------------------------LM_Iter %d --------------------------- \n',LM_Iter);
    fprintf('Normalized Param: %1.6f \n', norm_park0);
    
    % Run sim with param update
    Y_SIM_CELL = cell(num_inputs,1);
    park0 = save_param_org(:,LM_Iter); %norm_to_origin(norm_park0,B,selection_vector);
    fprintf('Actual Param for simulation : %1.6e \n', park0);
    
    SensFlag = 0;
    parfor idx = 1:num_inputs
        [Y_SIM_CELL{idx},~] = DFN_sim_casadi(p,exp_num{idx},Current_exp{idx}, Time_exp{idx}, Voltage_exp{idx}, T_amb{idx}, SensSelec, park0,SensFlag,Rc{idx});
    end
    
    v_sim = cell2mat(Y_SIM_CELL);
    
    % Save Info.
    y_minus_yfit = v_dat - v_sim;
    save_y_minus_yfit(:,LM_Iter) = y_minus_yfit;
    save_y_sim(:,LM_Iter) = v_sim;
    sigma_hat_sq = (1/(total_NT-num_param)) * (y_minus_yfit)'*(y_minus_yfit);
    
    % Save Various Metrics & Exit Criteria
    save_L2norm(LM_Iter) = norm(v_sim-v_dat,2);
    save_RMSE(LM_Iter) = rmse(v_dat,v_sim);
    chi_sq(LM_Iter) = (y_minus_yfit)'*W*(y_minus_yfit);
    redX2(LM_Iter) = chi_sq(LM_Iter) / (total_NT - num_param +1);
    
    save_param_exit(LM_Iter) = max(abs(delta_theta)); %[ZTG Change], removed /norm since delta_theta already norm
    save_chi_sq_exit(LM_Iter) = abs((chi_sq(LM_Iter) - chi_sq(LM_Iter-1)))/  chi_sq(LM_Iter-1); %[ZTG Change], checks whether cost function is converging
    save_chi_sq_AbsTol(LM_Iter) = abs(chi_sq(LM_Iter) - chi_sq(LM_Iter-1)); % Abs. Tol for Cost Function
    %%% Only simulation case.
    %     W = 1;
    %%% Experiment case.
    %     W = 1/(sigma_hat_sq);
    %%% Default version.
    %     W = 1/(y_dat'*y_dat);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %2. Check Exit condition
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% Multi-objective exit condition [ZTG change]
    % If (Cost Function decreases) AND (parameters converage OR
    % cost function converges)
    if (chi_sq(LM_Iter) < chi_sq(LM_Iter-1)) && ((save_chi_sq_AbsTol(LM_Iter) < chi_sq_Abs_exit_thresh) || (save_param_exit(LM_Iter) < param_exit_thresh) || (save_chi_sq_exit(LM_Iter) < chi_sq_exit_thresh))
        fprintf('Converged in parameters at %d iterations \n',LM_Iter)
        break;
    end
    
    % Exit conditions for max L-M iterations. Note: need if/else
    % because in some cases, on the last iteration the cost
    % function is actually worse. When this happens, need to throw
    % out most recent parameter set and revert to the best one.
    if LM_Iter == maxIter && (chi_sq(LM_Iter) > chi_sq(LM_Iter-1))
        fprintf('Max L-M iterations (%d) reached. Cost function fit worse. Throw out most recent parameters. \n',LM_Iter)
        
        % Reject current save_param(:,LM_Iter), back to previous one.
        park0 = save_param_org(:,LM_Iter-1);
        norm_park0 = save_param_nmz(:,LM_Iter-1);
        
        % Store best iteration for metric being tracked
        [rmse_best,ind_best] = min(paramID_out.save_RMSE); % determine best rmse value and corresponding ind.
        save_RMSE(LM_Iter) = rmse_best;
        
        chi_sq(LM_Iter) = chi_sq(ind_best);
        save_L2norm(LM_Iter) = save_L2norm(ind_best);
        save_delta_matrix(LM_Iter) = save_delta_matrix(ind_best);
        save_param_nmz(LM_Iter) = save_param_nmz(ind_best);
        save_param_org(LM_Iter) = save_param_org(ind_best);
        save_y_minus_yfit(LM_Iter) = save_y_minus_yfit(ind_best);
        save_y_sim(LM_Iter) = save_y_sim(ind_best);
        save_param_exit(LM_Iter) = save_param_exit(ind_best);
        save_chi_sq_exit(LM_Iter) = save_chi_sq_exit(ind_best);
        save_chi_sq_AbsTol(LM_Iter) = save_chi_sq_AbsTol(ind_best);
        save_lambda_matrix(LM_Iter) = save_lambda_matrix(ind_best);
        break;
    end
    
    % MaxIter but last iteration is acceptable
    if LM_Iter == maxIter && (chi_sq(LM_Iter) < chi_sq(LM_Iter-1))
        fprintf('Max L-M iterations (%d) reached. Cost function fit improved. \n',LM_Iter)
        break;
    end
    %         % Convergence in parameters
    %         %%%*******
    %         if(max(abs((delta_theta)./norm_park0)) < 1e-3)
    %             fprintf('Converged in parameters at %d iterations \n',LM_Iter)
    %             break;
    %         end
    %
    %          % Convergence in gradients
    %         if(max(abs((J_LM'*W*(y_dat-y_sim)))) < 1e-3)
    %             fprintf('Converged in gradients at %d iterations \n',LM_Iter)
    %             break;
    %         end
    %
    %         % Convergence in chi_sq
    %         if(save_L2norm(LM_Iter) < 1e-1)
    %             fprintf('Converged in reduced Chi-square at %d iterations \n',LM_Iter)
    %             break;
    %         end
    %
    %
    %         % Convergence in RMSE (<10mV)
    %         if(save_RMSE(LM_Iter) < 5e-3)
    %             fprintf('Converged in RMSE (less than 5mV) at %d iterations \n',LM_Iter)
    %             break;
    %         end
    %
    %         % Max iterations
    %         if LM_Iter == maxIter
    %             fprintf('Max iterations reached, %d iterations\n', LM_Iter);
    %             break;
    %         end
    
    % Display info.
    fprintf('Chi_sq: %e \n',chi_sq(LM_Iter));
    fprintf('Reduced chi_sq: %e \n',redX2(LM_Iter));
    fprintf('L2-norm: %f \n',save_L2norm(LM_Iter));
    fprintf('RMSE: %f \n',rmse(v_dat,v_sim));
    fprintf('Parameter convergence criterion: %f \n',save_param_exit(LM_Iter));%[ZTG Change]
    fprintf('Cost function rel. tolerance criterion: %f \n',save_chi_sq_exit(LM_Iter));%[ZTG Change]
    fprintf('Cost function abs. tolerance criterion: %f \n',save_chi_sq_AbsTol(LM_Iter));%[ZTG Change]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 3. Determine sensitivity calculation: Recalculate when cost
    %    function fit gets worse
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if(chi_sq(LM_Iter) > chi_sq(LM_Iter-1))
        
        % Reject current save_param(:,LM_Iter), back to previous one.
        park0 = save_param_org(:,LM_Iter-1);
        norm_park0 = save_param_nmz(:,LM_Iter-1);
        
        save_param_org(:,LM_Iter) = save_param_org(:,LM_Iter-1);
        save_param_nmz(:,LM_Iter) = save_param_nmz(:,LM_Iter-1);
        
        % Reject current chi_sq / RMSE as well. This way, cost function
        % will never get a worse fit than it begins with.
        chi_sq(LM_Iter) = chi_sq(LM_Iter-1);
        save_RMSE(LM_Iter) = save_RMSE(LM_Iter-1);
        
        % If cost function is worse in consecutive iterations, no
        % need to recalculate sensitivity because the parameters
        % have not changed. Only need to change lambda
        if SensCheck == 1
            fprintf('########## Recalculate Sensitivity at LM_Iter: %d ########## \n', LM_Iter)
            
            % Recalculate sensitivity
            V_CELL = cell(num_inputs,1);
            S_CELL = cell(num_inputs,1);
            
            SensFlag = 1;
            parfor idx = 1:num_inputs
                [V_CELL{idx}, ~, S_CELL{idx}] = DFN_sim_casadi(p,exp_num{idx},Current_exp{idx}, Time_exp{idx}, Voltage_exp{idx}, T_amb{idx}, SensSelec, park0, SensFlag,Rc{idx});
            end
            v_sim = cell2mat(V_CELL);
            normalized_sens_bar = origin_to_norm('sens',park0,bounds,selection_vector);
            S_LM = cell2mat(S_CELL);
            J_LM = bsxfun(@times,normalized_sens_bar', S_LM);
            
            clear V_LM_CELL
            clear S_LM_CELL
            
            % Recalculate lambda based on new sensitivity J_LM
            normalize_lambda = diag(J_LM'*W*J_LM);
            lambda = ctrl_lambda * diag(normalize_lambda);
            
            % Update parameter (INCREASE lambda to DECREASE parameter
            % update step size) -- Small steps in steepest-descent dir
            lambda = lambda*10;
            delta_theta = (lambda*bigI + J_LM'*W*J_LM)\(J_LM')*W*(v_dat - v_sim);
            save_delta_matrix(:,LM_Iter) = delta_theta;
            save_lambda_matrix{LM_Iter} = lambda;
            
            % Sensitivity was just recalculated; if cost function
            % worse in consecutive iteration(s), don't need to
            % recalculate sensitivity
            SensCheck = 0;
            
            LM_logic = vertcat(LM_logic,0);% variable to track what actions L-M is taking (-1: decrease lambda; 0: recalc sensitivity; 1: increase lambda)
            
        else
            fprintf('########## Increase Lambda at LM_Iter: %d ########## \n', LM_Iter)
            % Update parameter (INCREASE lambda to DECREASE parameter
            % update step size) -- Small steps in steepest-descent dir
            lambda = lambda*10;
            delta_theta = (lambda*bigI + J_LM'*W*J_LM)\(J_LM')*W*(v_dat - v_sim);
            save_delta_matrix(:,LM_Iter) = delta_theta;
            save_lambda_matrix{LM_Iter} = lambda;
            
            LM_logic = vertcat(LM_logic,1); % variable to track what actions L-M is taking (-1: decrease lambda; 0: recalc sensitivity; 1: increase lambda)
        end
    else %(chi_sq(LM_Iter) < chi_sq(LM_Iter-1))
        fprintf('=>=> Keep going this direction at LM_Iter: %d =>=> \n', LM_Iter)
        
        % Update parameter (DECREASE lambda to INCREASE parameter
        % update step size) -- Large steps (Gauss-Newton)
        lambda = lambda*0.1;
        delta_theta = (lambda*bigI + J_LM'*W*J_LM)\(J_LM')*W*(v_dat - v_sim);
        save_delta_matrix(:,LM_Iter) = delta_theta;
        save_lambda_matrix{LM_Iter} = lambda;
        SensCheck = 1; % turn on/keep on
        
        LM_logic = vertcat(LM_logic,-1);% variable to track what actions L-M is taking (-1: decrease lambda; 0: recalc sensitivity; 1: increase lambda)
    end
    
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


disp('After LM estimation')
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