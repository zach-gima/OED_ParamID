% DFN Voltage Obj. Fcn for parameter identification using fmincon

% By: ZTG 2019-5-3
% Purpose: Take in current parameter estimate and map to an RMSE value
% Inputs: (1) Initial parameter estimate
%         (2) Input profile

function [ cost,fgrad ] = V_obj(selected_params,selection_vector, Inputs, p)
   %% Load Input and Set Variables   
    Current_exp = Inputs.Current_exp;
    Time_exp = Inputs.Time_exp;
    Voltage_exp = Inputs.V_LM_CELL;
    T_amb = Inputs.T_amb_sim; % note comes in celcius
    exp_num = Inputs.exp_num;
    num_inputs = length(Voltage_exp);
    
    % In experimental ID, Rc needs to be identified for each experiment
    % (Rc_tune) and the Rc value must be attached to each profile
    if isfield(Inputs,'Rc')
        Rc = Inputs.Rc;
    else %M2M case
        Rc = cell(length(Current_exp),1); % Create a cell num_exp x 1
        Rc(:,1) = {p.R_c}; % For M2M case, just use nominal Rc value 
    end
    
    %% Simulate DFN
    V_LM_CELL = cell(num_inputs,1);
    S_LM_CELL = cell(num_inputs,1);
    
    % Simulate DFN and Calculate Sensitivties for initial parameter values
    SensFlag = 1;
    parfor idx = 1:num_inputs        
        [V_LM_CELL{idx}, ~, S_LM_CELL{idx}] = DFN_sim_casadi(p,Current_exp{idx}, Time_exp{idx}, Voltage_exp{idx}, T_amb{idx}, selection_vector, selected_params,SensFlag,Rc{idx});
    end
    
    %% Compute Cost
    v_sim = cell2mat(V_LM_CELL);
    sens = cell2mat(S_LM_CELL);
    v_dat = cell2mat(Voltage_exp); % truth/measured voltage data we're fitting to   

    N = length(v_dat);
    MSE = mean((v_dat - v_sim).^2);    
    RMSE = sqrt(MSE);

    fgrad_MSE = -2/N*sens'*(v_dat - v_sim); % MSE gradient
    fgrad = 1/(2*sqrt(MSE))*fgrad_MSE; %RMSE gradient

    %calculate the output measurement (rmse)
    cost = RMSE; %MSE
    fprintf('Voltage RMSE = %1.6f \n \n',cost);
end