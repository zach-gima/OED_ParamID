%% OED-CVX Parameter Estimation Function 
% Adapted from Saheong Park Exp_param_ID_normalized_saehong_G1G2G3G4.m
% By Zach Gima, 2018-2-22

function [park0, paramID_out, LM_Iter] = Param_ID(p,bounds,sel_k,selection_vector,theta_0,Inputs,filename_output,LM_options)
    %% For Savio

    % addpath('/global/home/users/sspark/casadi.bin')

    % parpool()

    % For Savio
    % all the cores on a node
    % parpool(str2num(getenv('SLURM_CPUS_ON_NODE'))); % 
    % many cores as requested.
    % parpool(str2num(getenv('SLURM_CPUS_PER_TASK')));

    %% Setup Parameter Variables
    % Parse L-M Conditions
    param_exit_thresh = LM_options.exit_cond(1);
    chi_sq_exit_thresh = LM_options.exit_cond(2);
    
    maxIter = LM_options.maxIter;
    ctrl_lambda = LM_options.ctrl_lambda;
    
    % Setup vectors related to parameter sensitivity selection
    SensSelec = selection_vector;
    Selected_params = theta_0(sel_k);  % Goes from 21x1 vector to 18x1 (in all params selected scenario)
    num_param = size(Selected_params,1); % number of param  

    
    LM_logic = 0; % variable to track what actions L-M is taking (-1: decrease lambda; 0: recalc sensitivity; 1: increase lambda)    

    % Check things are imported properly
    if ( size(nonzeros(SensSelec),1) ~= size(Selected_params,1)); error('size is wrong'); end

    %% Load Input and Set Variables
   
    Current_exp = Inputs.Current_exp;
    Time_exp = Inputs.Time_exp;
    Voltage_exp = Inputs.V_LM_CELL;
    T_amb = Inputs.T_amb_sim; % note comes in celcius
    
    num_inputs = length(Voltage_exp);
    y_dat = cell2mat(Voltage_exp); % Truth/experimental data
    total_NT  = size(y_dat,1);

    %%% *******************For Model-to-Model Comparison, don't need Rc value
    %%% Only important for matching IR drop between experimental data
    %%% and simulation. In which case, use Rc for each experiment
    % Rc = 0.0027; 

    %% Q (Covariance value)
    % [ZTG Change]
    % For Model-to-Model, set Q = 1
    Qfinal = 1;

    % Cell_Inputs = cell(num_inputs,1);
    % for i = 1: num_inputs
    %     Cell_Inputs{i} = eval(['Current_exp' num2str(i)]);
    % end
    % Concat_input = cell2mat(Cell_Inputs);
    % 
    % IntI = trapz(Concat_input);
    % Qchg=0.000062382391079*sqrt(trapz(Concat_input.^2))-0.004073969403817;
    % Qdchg=exp((0.017020092199588*sqrt(trapz(Concat_input.^2))))*exp(-8.466269267442632);
    % if IntI<-3500
    %     Qfactor=Qdchg;
    % elseif IntI>3500
    %     Qfactor=Qchg;
    % else
    %     Qfactor=Qchg/2+Qdchg/2;
    % end
    % 
    % Qfinal = Qfactor.^2;

    %% Run DFN or Sensitivity for Initial Parameter Guesses

    V_LM_CELL = cell(num_inputs,1);
    S_LM_CELL = cell(num_inputs,1);

    % [ZTG Change] Only need below if using experimental data (not
    % model-to-model)
    % Current_exp = cell(num_inputs,1);
    % Time_exp = cell(num_inputs,1);
    % Voltage_exp = cell(num_inputs,1);
    % Rc_exp = cell(num_inputs,1);

    % for i=1:num_inputs
    %     Current_exp{i} = eval(['Current_exp' num2str(i)]);
    %     Time_exp{i} = eval(['Time_exp' num2str(i)]);
    %     Voltage_exp{i} = eval(['Voltage_exp' num2str(i)]);
    %     Rc_exp{i} = eval(['Rc_exp' num2str(i)]);
    % end

    % Simulate DFN and Calculate Sensitivties for initial parameter values
    SensFlag = 1;
    parfor idx = 1:num_inputs
        [V_LM_CELL{idx}, ~, S_LM_CELL{idx}] = DFN_sim_casadi(p,Current_exp{idx}, Time_exp{idx}, Voltage_exp{idx}, T_amb{idx}, SensSelec, Selected_params,SensFlag);
    end

    V_LM = cell2mat(V_LM_CELL);
    S_LM = cell2mat(S_LM_CELL);

    clear V_LM_CELL
    clear S_LM_CELL

    % Calculate Information matrix
    fprintf('Determinant of information matrix: %e \n',det(S_LM'*S_LM));

    %% Normalization

    normalized_theta_bar = origin_to_norm('param',Selected_params,bounds,selection_vector);                     

%     Min_theta_bar = zeros(21,1);
%     Max_theta_bar = ones(21,1);
    Min_theta_bar = zeros(25,1);
    Max_theta_bar = ones(25,1);
    normalize_params_min = Min_theta_bar(sel_k);
    normalize_params_max = Max_theta_bar(sel_k);

    normalized_sens_bar = origin_to_norm('sens',Selected_params,bounds,selection_vector);
   
    %% Levenberg_Marquardt
    
    %create template output filename for saving each iteration of L-M data
    LM_filename_output = strcat(filename_output(1:end-4),'_iter_');

    % Sens -> LM variables
    y_sim = V_LM;
    %J_LM = (normalized_sens_bar)'.*S_LM; % Infeasible for Savio
    J_LM = bsxfun(@times,normalized_sens_bar',S_LM);
    fprintf('Determinant of information matrix(normalized): %e \n',det(J_LM'*J_LM));

    norm_park0 = normalized_theta_bar;
    origin_park0 = Selected_params;

    disp('---Before LM estimation---')
    fprintf('Normalized Param: %1.6f \n', norm_park0);
    fprintf('Actual Param: %1.6e \n', origin_park0);

    diff=norm(y_sim-y_dat,2);
    fprintf('Voltage difference L-2 norm before ID: %f \n',diff);

    % Determine initial RMSE [ZTG Note]
    y_minus_yfit = y_dat - y_sim;

    % W = 1/((y_dat)'*(y_dat)); % by default

    %%% Only simulation case.
    % W = 1

    %%% Experiment case.
    % sigma_hat_sq = (1/(total_NT-num_param)) * (y_minus_yfit)'*(y_minus_yfit);
    % W = 1/(sigma_hat_sq);

    %%% Default version.
    % W = 1/(y_dat'*y_dat);

    %%% Output error variance.
    W = 1/Qfinal;

    normalize_lambda = diag(J_LM'*W*J_LM);
    lambda = ctrl_lambda*diag(normalize_lambda); % When initial guess is far from true, set lambda BIG!

    bigI = eye(size(J_LM,2));
    delta_theta = ((lambda*bigI + J_LM'*W*J_LM)\(J_LM')*W*(y_dat - y_sim)); % Parameter Perturbation

    save_param_nmz = zeros(num_param,1);
    save_param_nmz(:,1) = norm_park0;
    save_param_org = zeros(num_param,1);
    save_param_org(:,1) = origin_park0;

    save_y_minus_yfit = zeros(total_NT,1);
    save_y_minus_yfit(:,1) = y_minus_yfit;
    save_y_sim(:,1) = y_sim;

    save_delta_matrix = zeros(num_param,1);
    save_delta_matrix(:,1) = delta_theta; %.* (normalize_sens');
    
    save_lambda_matrix{1} = lambda;
    
    % Calculate Exit Criteria for Initial Iteration
    chi_sq(1) = (y_minus_yfit)'*W*(y_minus_yfit);
    redX2(1) = chi_sq(1)/(total_NT - num_param +1);
    save_param_exit = []; %[ZTG Change]
    save_L2norm(1) = norm(y_sim-y_dat,2);
    save_RMSE(1) = rmse(y_dat,y_sim);
    save_chi_sq_exit = []; % for initial iteration, there isn't a previous iteration to calculate convergence with

    fprintf('Initial Reduced chi_sq: %e \n',redX2(1));
    fprintf('Initial L2-norm: %f \n',save_L2norm(1));
    fprintf('Initial RMSE: %f \n',save_RMSE(1));
    %     fprintf('Initial Parameter convergence criterion: %f \n',save_param_exit(1));%[ZTG Change]
    %     fprintf('Initial Cost function convergence criterion: %f \n',save_chi_sq_exit(1));%[ZTG Change]
    
    % If cost function is worse in consecutive iterations, no
    % need to recalculate sensitivity because the parameters
    % have not changed. Only need to change lambda
    SensCheck = 1; %initialize as on
    
    %For true initial parameter guess scenario, the exit conditions should
    %be met immediately; otherwise, L-M will push param values away before
    %checking exit conditions
    % For true initial params, RMSE ~1e-5, param_exit
    % ~ 2e-3
    if(max(abs(delta_theta)) < 1e-4) && (save_RMSE(1) <1e-4)%(save_chi_sq_exit(1) < 1e-3)
        LM_Iter = 1;
        park0 = save_param_org(:,LM_Iter); %norm_to_origin(norm_park0,B,selection_vector);
        fprintf('Converged in parameters at %d iterations \n',LM_Iter)
        fprintf('Actual Param for simulation : %1.6e \n', park0);
    else
        for LM_Iter=2:maxIter

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % 1. Param updates & Simulation Run
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

            norm_park0 = norm_park0 + delta_theta;
            norm_park0 = min(max(normalize_params_min,norm_park0),normalize_params_max);
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
                [Y_SIM_CELL{idx},~] = DFN_sim_casadi(p,Current_exp{idx}, Time_exp{idx}, Voltage_exp{idx}, T_amb{idx}, SensSelec, park0,SensFlag);
            end

            y_sim = cell2mat(Y_SIM_CELL);

            % Save Info.
            y_minus_yfit = y_dat - y_sim;
            save_y_minus_yfit(:,LM_Iter) = y_minus_yfit;
            save_y_sim(:,LM_Iter) = y_sim;
            sigma_hat_sq = (1/(total_NT-num_param)) * (y_minus_yfit)'*(y_minus_yfit);

            % Save Various Metrics & Exit Criteria
            save_L2norm(LM_Iter) = norm(y_sim-y_dat,2);
            save_RMSE(LM_Iter) = rmse(y_dat,y_sim);
            chi_sq(LM_Iter) = (y_minus_yfit)'*W*(y_minus_yfit);
            redX2(LM_Iter) = chi_sq(LM_Iter) / (total_NT - num_param +1);

            save_param_exit(LM_Iter) = max(abs(delta_theta)); %[ZTG Change], removed /norm since delta_theta already norm
            save_chi_sq_exit(LM_Iter) = abs((chi_sq(LM_Iter) - chi_sq(LM_Iter-1)))/  chi_sq(LM_Iter-1); %[ZTG Change], checks whether cost function is converging
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
            % If change in parameters converage and change in error converges 
            if(save_param_exit(LM_Iter) < param_exit_thresh) || (save_chi_sq_exit(LM_Iter) < chi_sq_exit_thresh)
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
             fprintf('RMSE: %f \n',rmse(y_dat,y_sim));
             fprintf('Parameter convergence criterion: %f \n',save_param_exit(LM_Iter));%[ZTG Change]
             fprintf('Cost function convergence criterion: %f \n',save_chi_sq_exit(LM_Iter));%[ZTG Change]

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
                    V_LM_CELL = cell(num_inputs,1);
                    S_LM_CELL = cell(num_inputs,1);

                    SensFlag = 1;
                    parfor idx = 1:num_inputs
                        [V_LM_CELL{idx}, ~, S_LM_CELL{idx}] = DFN_sim_casadi(p,Current_exp{idx}, Time_exp{idx}, Voltage_exp{idx}, T_amb{idx}, SensSelec, park0, SensFlag);
                    end
                    y_sim = cell2mat(V_LM_CELL);
                    normalized_sens_bar = origin_to_norm('sens',park0,bounds,selection_vector);
                    S_LM = cell2mat(S_LM_CELL);
                    J_LM = bsxfun(@times,normalized_sens_bar', S_LM);

                    clear V_LM_CELL
                    clear S_LM_CELL
                    
                    % Recalculate lambda based on new sensitivity J_LM
                    normalize_lambda = diag(J_LM'*W*J_LM);
                    lambda = ctrl_lambda * diag(normalize_lambda);
                    
                    % Update parameter (INCREASE lambda to DECREASE parameter
                    % update step size) -- Small steps in steepest-descent dir
                    lambda = lambda*10;
                    delta_theta = (lambda*bigI + J_LM'*W*J_LM)\(J_LM')*W*(y_dat - y_sim);
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
                    delta_theta = (lambda*bigI + J_LM'*W*J_LM)\(J_LM')*W*(y_dat - y_sim);
                    save_delta_matrix(:,LM_Iter) = delta_theta;
                    save_lambda_matrix{LM_Iter} = lambda;
                    
                    LM_logic = vertcat(LM_logic,1); % variable to track what actions L-M is taking (-1: decrease lambda; 0: recalc sensitivity; 1: increase lambda)
                end
            else %(chi_sq(LM_Iter) < chi_sq(LM_Iter-1))
                fprintf('=>=> Keep going this direction at LM_Iter: %d =>=> \n', LM_Iter)
                
                % Update parameter (DECREASE lambda to INCREASE parameter
                % update step size) -- Large steps (Gauss-Newton)
                lambda = lambda*0.1;
                delta_theta = (lambda*bigI + J_LM'*W*J_LM)\(J_LM')*W*(y_dat - y_sim);
                save_delta_matrix(:,LM_Iter) = delta_theta;
                save_lambda_matrix{LM_Iter} = lambda;
                SensCheck = 1; % turn on/keep on 
                
                LM_logic = vertcat(LM_logic,-1);% variable to track what actions L-M is taking (-1: decrease lambda; 0: recalc sensitivity; 1: increase lambda)
            end
                        
            % Save ParamID results every iteration
            paramID_out.Time_exp = Time_exp;
            paramID_out.Current_exp = Current_exp;
            paramID_out.Voltage_exp = Voltage_exp;
            paramID_out.save_chi_sq = chi_sq;
            paramID_out.save_L2norm = save_L2norm;
            paramID_out.save_RMSE = save_RMSE;
            paramID_out.save_delta_matrix = save_delta_matrix;
            paramID_out.save_param_nmz = save_param_nmz;
            paramID_out.save_param_org = save_param_org;
            paramID_out.save_y_minus_yfit = save_y_minus_yfit;
            paramID_out.save_y_sim = save_y_sim;
            paramID_out.save_param_exit = save_param_exit;
            paramID_out.save_chi_sq_exit = save_chi_sq_exit;
            paramID_out.save_lambda_matrix = save_lambda_matrix;
            paramID_out.y_dat = y_dat;    
            paramID_out.LM_logic = LM_logic;
            
            save(strcat(LM_filename_output,num2str(LM_Iter),'.mat'),'park0','paramID_out','LM_Iter');
        end
    end

    disp('After LM estimation')
    fprintf('Final Parameter Values: %1.6f \n',park0)    
    
    %% Concatenate Outputs & Save Results
    paramID_out.Time_exp = Time_exp;
    paramID_out.Current_exp = Current_exp;
    paramID_out.Voltage_exp = Voltage_exp;
    paramID_out.save_chi_sq = chi_sq;
    paramID_out.save_L2norm = save_L2norm;
    paramID_out.save_RMSE = save_RMSE;
    paramID_out.save_delta_matrix = save_delta_matrix;
    paramID_out.save_param_nmz = save_param_nmz;
    paramID_out.save_param_org = save_param_org;
    paramID_out.save_y_minus_yfit = save_y_minus_yfit;
    paramID_out.save_y_sim = save_y_sim;
    paramID_out.save_param_exit = save_param_exit;
    paramID_out.save_chi_sq_exit = save_chi_sq_exit;
    paramID_out.save_lambda_matrix = save_lambda_matrix;
    paramID_out.y_dat = y_dat;
    paramID_out.LM_logic = LM_logic;

end