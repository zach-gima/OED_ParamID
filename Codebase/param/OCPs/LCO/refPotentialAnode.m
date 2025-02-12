%% Reference Potential for Neg Electrode: Unref(theta_n)
%   Created July 12, 2011 by Scott Moura

function [Uref,varargout] = refPotentialAnode(p,theta)

% Commented for CasADi

% % if(~isreal(theta))
% %     beep;
% %     error('dfn:err','Complex theta_n');
% % %     pause;
% % end

% % Polynomail Fit
% Uref = ppvalFast(p.Uppn,theta);

% DUALFOIL: MCMB 2528 graphite (Bellcore) 0.01 < x < 0.9
Uref = 0.194+1.5*exp(-120.0*theta) ...
     +0.0351*tanh((theta-0.286)/0.083) ... 
     - 0.0045*tanh((theta-0.849)/0.119) ...
     - 0.035*tanh((theta-0.9233)/0.05) ...
     - 0.0147*tanh((theta-0.5)/0.034) ...
     - 0.102*tanh((theta-0.194)/0.142) ...
     - 0.022*tanh((theta-0.9)/0.0164) ...
     - 0.011*tanh((theta-0.124)/0.0226) ...
     + 0.0155*tanh((theta-0.105)/0.029);

% Gradient of OCP wrt theta
if(nargout >= 2)

%     % Polynomial Fit
%     dUref = ppvalFast(p.dUppn,theta);
%     varargout{1} = dUref / p.c_s_n_max;

dUref = -1.5*(120.0/p.c_s_n_max)*exp(-120.0*theta)  ...
 +(0.0351/(0.083*p.c_s_n_max))*((cosh((theta-0.286)/0.083)).^(-2)) ...
 -(0.0045/(p.c_s_n_max*0.119))*((cosh((theta-0.849)/0.119)).^(-2)) ...
 -(0.035/(p.c_s_n_max*0.05))*((cosh((theta-0.9233)/0.05)).^(-2)) ...
 -(0.0147/(p.c_s_n_max*0.034))*((cosh((theta-0.5)/0.034)).^(-2)) ...
 -(0.102/(p.c_s_n_max*0.142))*((cosh((theta-0.194)/0.142)).^(-2)) ...
 -(0.022/(p.c_s_n_max*0.0164))*((cosh((theta-0.9)/0.0164)).^(-2)) ...
 -(0.011/(p.c_s_n_max*0.0226))*((cosh((theta-0.124)/0.0226)).^(-2)) ...
 +(0.0155/(p.c_s_n_max*0.029))*((cosh((theta-0.105)/0.029)).^(-2));
varargout{1} = dUref;

end

% Gradient of OCP wrt temperature
if(nargout >= 3)
    
    dUdT = 0;
    varargout{2} = dUdT;
    
end

