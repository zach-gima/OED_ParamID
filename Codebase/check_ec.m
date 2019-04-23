% function for checking parameter identification routine exit conditions
function [exit_logic] = check_ec(v_dat,v_sim,delta_theta,Iter,SCD_options)

    % Parse Optimization Routine options
    maxIter = SCD_options.maxIter;

    % SCD_options.exit_cond = [param_exit_thresh, chi_sq_rel_thresh, chi_sq_abs_thresh];
    param_exit_thresh = SCD_options.param_exit_thresh;
    chi_sq_rel_thresh = SCD_options.chi_sq_rel_thresh;
    chi_sq_abs_thresh = SCD_options.chi_sq_abs_thresh;
    % Compute exit condition-relevant metrics
    y_minus_yfit = v_dat - v_sim;
    
    % Save Various Metrics & Exit Criteria
%     save_RMSE = rmse(v_dat,v_sim);
    chi_sq = (y_minus_yfit)'*W*(y_minus_yfit);
    
    param_exit = max(abs(delta_theta)); % Convergence in the parameter estimates
    chi_sq_RelTol = abs((chi_sq - chi_sq(Iter-1)))/  chi_sq(Iter-1); % Rel. Tol for Cost Function
    chi_sq_AbsTol = abs(chi_sq - chi_sq(Iter-1)); % Abs. Tol for Cost Function
 
    % Display info.
    fprintf('Chi_sq: %e \n',chi_sq);
    fprintf('RMSE: %f \n',rmse(v_dat,v_sim));
    fprintf('Parameter convergence criterion: %f \n',param_exit);%[ZTG Change]
    fprintf('Cost function rel. tolerance criterion: %f \n',chi_sq_RelTol);%[ZTG Change]
    fprintf('Cost function abs. tolerance criterion: %f \n',chi_sq_AbsTol);%[ZTG Change]
    
    %%% Multi-objective exit condition [ZTG change]
    % If (Cost Function decreases) AND (parameters converage OR
    % cost function converges)
    if (chi_sq < chi_sq(Iter-1)) && ((chi_sq_AbsTol < chi_sq_abs_thresh) || (param_exit < param_exit_thresh) || (chi_sq_RelTol < chi_sq_rel_thresh))
        fprintf('Converged in parameters at %d iterations \n',Iter)
        exit_logic = true;
    end
    
    if Iter == maxIter
       frpintf('Max Iterations Reached')
       exit_logic = true;
    end
    
    exit_logic = false;
  
end