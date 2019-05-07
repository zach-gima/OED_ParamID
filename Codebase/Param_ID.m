%% OED-CVX Parameter Estimation Function 
% Adapted from Saheong Park Exp_param_ID_normalized_saehong_G1G2G3G4.m
% By Zach Gima, 2018-2-22

function [theta_ID, paramID_out] = Param_ID(p,bounds,sel_k,selection_vector,theta_0,Inputs)
    %% Setup Parameter Variables
    
    % Setup vectors related to parameter sensitivity selection
    selected_params = theta_0(sel_k);  % Goes from 21x1 vector to 18x1 (in all params selected scenario)
    num_param = size(selected_params,1); % number of param  
    
    LM_logic = 0; % variable to track what actions L-M is taking (-1: decrease lambda; 0: recalc sensitivity; 1: increase lambda)    

    % Check things are imported properly
    if ( size(nonzeros(selection_vector),1) ~= size(selected_params,1)); error('size is wrong'); end

    
    %% fmincon
    lb = bounds.min(sel_k);
    ub = bounds.max(sel_k);
    
    opt = optimoptions('fmincon','Display','iter','Algorithm','sqp',...
    'SpecifyObjectiveGradient',true,'Diagnostics','on');  %'CheckGradients',true,'FiniteDifferenceType','central',
    
    % create anonymous function to pass in other parameters needed for V_obj
    fh = @(x)V_obj(x,batch_data,theta_sx_temp,p);
    
    [theta_ID,FVAL,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN]= ...
       fmincon(fh,theta.guess(ID_p.ID_idx),[],[],[],[],lb,ub,[],opt);
    
    
    %% Concatenate Outputs & Save Results
%     paramID_out.exp_num = exp_num;
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

end