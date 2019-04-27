function [p] = update_p(p,theta)
%Creadted By Dylan kato
%4/26/19
%Takes a vector of current parameters and inserts them into p

str={'p.D_s_n0', 'p.D_s_p0', 'p.R_s_n', 'p.R_s_p', 'p.epsilon_s_n', 'p.epsilon_s_p', 'p.sig_n', 'p.sig_p' ...
         'p.ElecFactorD', 'p.epsilon_e_n', 'p.epsilon_e_s', 'p.epsilon_e_p', 'p.ElecFactorK', 'p.t_plus', ...
         'p.ElecFactorDA','p.k_n0','p.k_p0', 'p.R_f_n', 'p.R_f_p','p.n_Li_s', 'p.c_e0', 'p.E.Dsn','p.E.Dsp','p.E.kn','p.E.kp'};
iter = [1:4,7:19,21:25]; %non eq structure params
for i = iter
    if str{i}(3)== 'E'
    p.E.(str{i}(5:end)) = theta(i);
    else
    p.(str{i}(3:end)) = theta(i);
    end
end
end