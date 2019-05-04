% function for checking parameter identification routine exit conditions
function [ec,exit_logic] = check_ec(v_dat,v_sim,ec,delta_theta,Iter,SCD_options,idx)

    % Parse Optimization Routine options
    maxIter = SCD_options.maxIter;
    param_exit_thresh = SCD_options.param_exit_thresh;
    chi_sq_rel_thresh = SCD_options.chi_sq_rel_thresh;
    chi_sq_abs_thresh = SCD_options.chi_sq_abs_thresh;
    
    % Compute exit condition-relevant metrics
    y_minus_yfit = v_dat - v_sim;
    
    % Save Various Metrics & Exit Criteria
    W = 1;
    ec.chi_sq(Iter) = (y_minus_yfit)'*W*(y_minus_yfit); % Chi squared for this iteration and parameter
    ec.chi_sq_mem(idx,Iter) = ec.chi_sq(Iter);  % Chi squared memory holds the prev. chi_sq value achieved for every parameter direction searched
    
    ec.param_exit(Iter) = max(abs(delta_theta)); % Convergence in the parameter estimates
    ec.chi_sq_RelTol(Iter) = max( abs((ec.chi_sq_mem(:,Iter) - ec.chi_sq_mem(:,Iter-1))) ./ ec.chi_sq_mem(:,Iter-1) ); % Rel. Tol for Cost Function
    ec.chi_sq_AbsTol(Iter) = abs(ec.chi_sq(Iter) - ec.chi_sq(Iter-1)); % Abs. Tol for Cost Function    
    
    % Display info.
    fprintf('Chi_sq: %e \n',ec.chi_sq(Iter));
    fprintf('RMSE: %f \n',rmse(v_dat,v_sim));
    fprintf('Parameter convergence criterion: %f \n',ec.param_exit(Iter));
    fprintf('Cost function rel. tolerance criterion: %f \n',ec.chi_sq_RelTol(Iter));
    fprintf('Cost function abs. tolerance criterion: %f \n',ec.chi_sq_AbsTol(Iter));
    
    % Evaluate multi-objective exit condition
    % If (Cost Function decreases) AND (parameters converage OR
    % cost function rel/abs converges)
    if (ec.chi_sq(Iter) < ec.chi_sq(Iter-1)) && ( (ec.param_exit(Iter) < param_exit_thresh) || (ec.chi_sq_AbsTol(Iter) < chi_sq_abs_thresh) || (ec.chi_sq_RelTol(Iter) < chi_sq_rel_thresh) )
        fprintf('Converged in parameters at %d iterations \n',Iter)
        exit_logic = true;
    end
    
    if Iter == maxIter
       fprintf('Max Iterations Reached')      
       exit_logic = true;
    end
    
    exit_logic = false;
  
end