%% Error analysis 
function [J_LM,ci95,sigma_y,covar_p,alg_states] = conf_interval(p, Inputs, SensSelec, park0)
    
    % Parse Inputs
    Current_exp = Inputs.Current_exp;
    Time_exp = Inputs.Time_exp;
    Voltage_exp = Inputs.V_LM_CELL;
    T_amb = Inputs.T_amb_sim;
    exp_num = Inputs.exp_num;

    num_inputs = length(Voltage_exp);
    y_dat = cell2mat(Voltage_exp); % Truth/experimental data
    total_NT  = size(y_dat,1);
    
    % In experimental ID, Rc needs to be identified for each experiment
    % (Rc_tune) and the Rc value must be attached to each profile
    if isfield(Inputs,'Rc')
        Rc = Inputs.Rc;
    else %M2M case
        Rc = cell(length(Current_exp),1);
        Rc(:,1) = {p.R_c}; % For M2M case, just use nominal Rc value 
    end
    
    % Calculate sensitivity at param_opt
    V_LM_CELL = cell(num_inputs,1);
    S_LM_CELL = cell(num_inputs,1);

    % Recalculate sensitivity (i.e. Jacobians) at final parameter values.
    % Necessary sometimes if the exit conditions are not met and the paramID
    % routine has to be manually stopped. In that case, the most recently
    % calculated Jacobian may not correspond to the identified parameter set.
    SensFlag = 1;
    parfor idx = 1:num_inputs
        [V_LM_CELL{idx}, alg_states{idx}, S_LM_CELL{idx}] = DFN_sim_casadi(p,exp_num{idx},Current_exp{idx}, Time_exp{idx}, Voltage_exp{idx},T_amb{idx}, SensSelec, park0, SensFlag, Rc{idx});
    end

    y_sim = cell2mat(V_LM_CELL);
    J_LM = cell2mat(S_LM_CELL);
    
%     y_minus_yfit = y_dat - y_sim;
%     rmse_yfit = rmse(y_dat,y_sim);

    % Normalizing J_LM before calculating conf. interval % [ZTG Change]
%     normalized_sens_bar = originS_to_normS(park0,bounds,selection_vector);
%     S_LM = cell2mat(S_LM_CELL);
%     J_LM_norm = bsxfun(@times,normalized_sens_bar', J_LM);
%     covar_p_norm = inv(J_LM'* weight * J_LM);
%     sigma_p_norm = sqrt(diag(covar_p));
    
    %%% Only simulation case.
    weight = 1; 

    %%% Experiment case.
%     DoF = total_NT - num_param +1;
%     sig_y_hat = (1/DoF)*(y_minus_yfit'*y_minus_yfit);
%     weight = eye(1)/sig_y_hat;

    %%% Q Case  
%     % Experimental Case
%     const = 2.903580414003080e-05;
%     Qfinal = zeros(num_inputs,1);
%     
%     W = cell(num_inputs,1);
%     for ii = 1:num_inputs
%         stdDev = const*norm(Current_exp{ii},2);
%         Qfinal(ii) = stdDev ^ 2;
%         W_temp = 1/Qfinal(ii);
%         W{ii} = ones(length(Current_exp{ii}),length(Current_exp{ii}))*W_temp;
%     end
%     W = sparse(blkdiag(W{:}));

    covar_p = inv(J_LM'* weight * J_LM);
    sigma_p = sqrt(diag(covar_p));

    % sigma_y = zeros(total_NT,1);
    % for i=1:NT
    %  sigma_y(i) = J_LM(i,:) * covar_p * J_LM(i,:)';        
    % end
    % sigma_y = sqrt(sigma_y);
    sigma_y = sqrt(diag(J_LM * covar_p * J_LM'));

    ci = 0.95;
    alpha = 1 - ci;
    T_multiplier = tinv(1-alpha/2, total_NT-1);
    
    % the multiplier is large here because there is so little data. That means
    % we do not have a lot of confidence on the true average or standard
    % deviation
    ci95 = T_multiplier*sigma_p/sqrt(total_NT);

    % disp('Confidence interval')
    % [park0 - ci95, park0 + ci95]

end