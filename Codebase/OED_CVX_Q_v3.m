clc;
clear all;
close all;
%% Optimal Experiment Design (OED) via Convex Programming

% Created Mar. 30. 2017 by Saehong Park

% Dependency: cvx
tic

fs = 18;
%%

%   UNCERTAIN PARAMETERS, theta
%   1  : D_s_n       => p.D_s_n0
%   2  : D_s_p       => p.D_s_p0
%   3  : R_s_n       => p.R_s_n
%   4  : R_s_p       => p.R_s_p
%   5  : sig_n       => p.sig_n
%   6  : sig_p       => p.sig_p
%   7  : D_e         => p.ElecFactorD --> multiply factor
%   8  : epsilon_e_n => p.epsilon_e_n
%   9  : epsilon_e_s => p.epsilon_e_s
%   10 : epsilon_e_p => p.epsilon_e_p
%   11 : kappa       => p.ElecFactorK --> multiply factor
%   12 : t_plus      => p.t_plus
%   13 : dactivity   => p.ElecFactorDA --> multiply factor
%   14 : k_n0        => p.k_n0
%   15 : k_p0        => p.k_p0
%   16 : R_f_n       => p.R_f_n
%   17 : R_f_p       => p.R_f_p
%   18 : c_e0        => p.c_e


params = {'$D_s^{^{\_}}$','$D_s^+$','$R_s^-$','$R_s^+$',...
    '$\sigma^{^{\_}}$','$\sigma^{+}$','$D_e$','$\varepsilon_e^{^{\_}}$','$\varepsilon_e^{sep}$','$\varepsilon_e^+$',...
    '$\kappa$','$t_{c}^{0}$','$\frac{d \ln f_{c/a}}{d \ln c_e}$','$k^{^{\_}}$','$k^+$','$R_f^{^{\_}}$','$R_f^+$','$c_{e_0}$'};

run ../param/params_NCA

% selvec = [1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1];
% selvec = [1;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0];
%%% DIVISION A
% selvec = [1;1;1;1;0;0;0;0;0;0;0;0;0;0;0;0;0;0]; % Tier 1
% selvec = [0;0;0;0;0;0;1;1;0;0;1;0;0;0;0;0;0;0]; % Tier 2
% selvec = [0;0;0;0;0;0;0;0;1;1;0;0;1;0;0;1;1;1]; % Tier 3
% selvec = [0;0;0;0;1;1;0;0;0;0;0;1;0;1;1;0;0;0]; % Tier 4

%%% DIVISION B (Since June 9)
% selvec = [0;0;0;1;0;0;1;0;0;0;1;0;0;0;0;0;0;0]; % Tier 1
% selvec = [0;0;1;0;0;0;0;1;0;0;0;0;0;0;0;0;0;0]; % Tier 2
% selvec = [1;1;0;0;0;0;0;0;1;0;0;0;0;0;0;0;0;0]; % Tier 3
% selvec = [0;0;0;0;0;0;0;0;0;1;0;0;1;0;0;1;1;1]; % Tier 4
% selvec = [0;0;0;0;1;1;0;0;0;0;0;0;0;1;1;0;0;0]; % Tier 5
% selvec = [0;0;0;0;0;0;0;0;0;0;0;1;0;0;0;0;0;0]; % Tier 6

%%% DIVISION C (Since June 15)
% selvec = [0;0;1;1;0;0;1;1;0;0;1;0;0;0;0;0;0;0]; % Tier 1
% selvec = [1;1;0;0;0;0;0;0;0;1;0;0;0;0;0;0;0;0]; % Tier 2
% selvec = [0;0;0;0;0;0;0;0;1;0;0;0;1;0;0;1;1;1]; % Tier 3
% selvec = [0;0;0;0;1;1;0;0;0;0;0;1;0;1;1;0;0;0]; % Tier 4
% selvec = [0;0;0;0;0;0;0;0;0;0;0;1;0;0;0;0;0;0]; % Tier 5

%%% DIVISION D (Since July 4)
% selvec = [0;0;1;1;0;0;0;1;0;1;0;0;0;0;0;0;0;0]; % Tier 1
% selvec = [1;1;0;0;0;0;1;0;0;0;1;0;0;0;0;0;0;0]; % Tier 2
% selvec = [0;0;0;0;0;0;0;0;1;0;0;0;1;0;0;1;1;1]; % Tier 3
% selvec = [0;0;0;0;1;1;0;0;0;0;0;1;0;1;1;0;0;0]; % Tier 4


%%% DIVISION E (Since Nov 2)
% selvec = [0;0;1;1;0;0;0;0;0;0;0;0;0;0;0;0;0;0]; % Tier 1
% selvec = [1;1;0;0;0;0;1;1;0;0;1;0;1;0;0;0;0;0]; % Tier 2
% selvec = [0;0;0;0;0;0;0;0;0;1;0;0;0;1;0;1;1;1]; % Tier 3
selvec = [0;0;0;0;1;1;0;0;1;0;0;1;0;0;1;0;0;0]; % Tier 4
% selvec = [1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1]; % ALL

%   UNCERTAIN PARAMETERS, theta
%   1  : D_s_n       => p.D_s_n0
%   2  : D_s_p       => p.D_s_p0
%   3  : R_s_n       => p.R_s_n
%   4  : R_s_p       => p.R_s_p
%   5  : sig_n       => p.sig_n
%   6  : sig_p       => p.sig_p
%   7  : D_e         => p.ElecFactorD --> multiply factor
%   8  : epsilon_e_n => p.epsilon_e_n
%   9  : epsilon_e_s => p.epsilon_e_s
%   10 : epsilon_e_p => p.epsilon_e_p
%   11 : kappa       => p.ElecFactorK --> multiply factor
%   12 : t_plus      => p.t_plus
%   13 : dactivity   => p.ElecFactorDA --> multiply factor
%   14 : k_n0        => p.k_n0
%   15 : k_p0        => p.k_p0
%   16 : R_f_n       => p.R_f_n
%   17 : R_f_p       => p.R_f_p
%   18 : c_e0        => p.c_e




sel_k = find(selvec);

%% Normalization

% % % % Log normalization
% % % % Multiply log(\theta_max / \theta_min) * theta_0

% % % % minmax normalization
% % % % Multiply (\theta_max - \theta_min)

% % % % [0 \theta_0]
% % % % Multiply norminal param

% % % % Parameter Min/Max values

Dsn_max = 1.05E-12; % LCO, since NCA literature is too small
Dsn_min = 2.25E-16; % LCO, ""
Dsp_max = 1.00E-12; % LCO, ""
Dsp_min = 2.00E-16; % LCO, ""

Rsn_min = 1.0e-06;
Rsn_max = 100e-06;
Rsp_min = 1.0e-06;
Rsp_max = 100e-06;

Sig_n_max = 500; % NCA, times 1 by nominal param
Sig_n_min = 50;  % NCA, times 1/2 by nominal param
Sig_p_max = 500; % NCA, times 5 by nominal param
Sig_p_min = 50;  % NCA, times 1/2 by nominal param

ElecFactorD_max = 1.5; % function, times 10 by nominal param
ElecFactorD_min = 0.5; % function, times 1/10 by nominal param

eps_e_n_max = 0.45;
eps_e_n_min = 0.18;
eps_e_s_max = 0.5;
eps_e_s_min = 0.45;
eps_e_p_max = 0.33;
eps_e_p_min = 0.18;

ElecFactorK_max = 1.5;
ElecFactorK_min = 0.5;

t_plus_max = 0.363;
t_plus_min = 0.36;

ElecFactorDA_max = 1.5;
ElecFactorDA_min = 0.5;

k_n0_max = 10*p.k_n0;
k_n0_min = 0.1*p.k_n0;
k_p0_max = 10*p.k_p0;
k_p0_min = 0.1*p.k_p0;

Rfn0_min = 1e-5;
Rfn0_max = 1e-3;
Rfp0_min = 1e-4;
Rfp0_max = 1e-3;

c_e0_max = 1500;
c_e0_min = 500;

normalize_param = [ (1/log10(exp(1)))*log10(Dsn_max / Dsn_min)*p.D_s_n0; ...
                    (1/log10(exp(1)))*log10(Dsp_max / Dsp_min)*p.D_s_p0; ...
                    Rsn_max - Rsn_min;
                    Rsp_max - Rsp_min;
                    (1/log10(exp(1)))*log10(Sig_n_max/Sig_n_min) * p.sig_n; ...
                    (1/log10(exp(1)))*log10(Sig_p_max/Sig_p_min) * p.sig_p; ...
                    (ElecFactorD_max - ElecFactorD_min); ...
                    eps_e_n_max - eps_e_n_min ; ...
                    eps_e_s_max - eps_e_s_min ; ...
                    eps_e_p_max - eps_e_p_min ; ...
                    (ElecFactorK_max - ElecFactorK_min) ; ...
                    (t_plus_max - t_plus_min) ; ...
                    (ElecFactorDA_max - ElecFactorDA_min); ...
                    (1/log10(exp(1)))*log10(k_n0_max / k_n0_min) * p.k_n0; ...
                    (1/log10(exp(1)))*log10(k_p0_max / k_p0_min) * p.k_p0; ...
                    Rfn0_max - Rfn0_min; ...
                    Rfp0_max - Rfp0_min; ...
                    c_e0_max - c_e0_min; ...
                    ];
normalize_param = normalize_param' ; % row maxtrix

normalize_param = normalize_param(sel_k);

%% Merging sensitivity result
%%% Saehong's workspace

% NumFiles = 18;
% 
% 
% file_path = 'SensResults/G3/';
% for i=1:NumFiles
%     data_path = [file_path 'sensG3_' num2str(i) '.mat'];
%     load(data_path); % stuct file: sens
%     savefile = 'sens';
%     
%     new_name = ['sensALL_' num2str(i+180+180+180+60+60+60) '.mat'];
%     save(new_name ,savefile)
% 
% end

%% Load Input Library

C = {}; % Row vector
cost = {}; % Column vector
Profile = {}; % DFN profile
QQQ={};
idx = 1; % index for pruned inputs.
index =[];
num_inputs = 738; % Total number of inputs


SOCint=zeros(num_inputs,1);
SOCend=zeros(num_inputs,1);


fn= 'SensResults/ALL/';
% fn= 'SensResults/G11new/';
disp(['Loaded Sensitivity directory:  ' fn]); 

for i=1:num_inputs
    data_path = [fn 'sensALL_' num2str(i) '.mat'];
    try
        load(data_path); % stuct file: sens
    catch
        fprintf(['Missing: ' data_path '\n']);
        sens.S3=1;
        sens.DFNout = 1;
        sens.Cost = 1e6;
        sens.Err_flag = 1;
    end
    
    S3 = sens.S3;
    if S3 ~= 1
        S3 = S3(:,sel_k); % Some input that has Err_flag=1, doesn't have matrix form. Ignore that
    end
    S3 = normalize_param .* S3;
    
    DFN_out = sens.DFNout;
    Cost = sens.COST;
    t= linspace(1,Cost,Cost);
    Nt = length(find(selvec));%size(S3,2);
    
    Err_flag = sens.Err_flag;
%     determinant=det(STS);
    
    if Err_flag == 0
        if max(DFN_out.T2) >= 273+35
            fprintf('Input %d violates temperature condition (T2) \n', i);
            Err_flag = 1;
        end
    end

    
    if Err_flag == 0
        %fprintf('%d\n',i);
        cost{idx} = Cost;
        Profile{idx} = DFN_out;
        Profile{idx}.index = i; % Save original input number
        
        % InitSOC, EndSOC
        t=1:length(DFN_out.Current);
        %fprintf('delta_soc for input %d: %f \n',i,trapz(t,DFN_out.Current*0.0744)/(3600*2.9));
        deltaSOC=-trapz(t,DFN_out.Current*0.0744)/(3600*2.9);
        V0=sens.DFNout.Volt(1);
        SOCint(idx)=and(le(V0,3.99),ge(V0,3.85))*.79300+and(le(V0,3.84999),ge(V0,3.68))*.592291+and(le(V0,3.6799),ge(V0,3.52))*.379585+...
        and(le(V0,3.51999),ge(V0,3.2))*.184100;
        SOCend(idx)=SOCint(idx)+deltaSOC;
        IntI=trapz(sens.DFNout.Current);
        
        Qchg=0.000062382391079*sqrt(trapz(sens.DFNout.Current.^2))-0.004073969403817;
       
        Qdchg=exp((0.017020092199588*sqrt(trapz(sens.DFNout.Current.^2))))*exp(-8.466269267442632);
        if IntI<-3500
            Qfactor=Qdchg;
        elseif IntI>3500
            Qfactor=Qchg;
        else
            Qfactor=Qchg/2+Qdchg/2;
        end
        QQQ{idx}=Qfactor;
        Q=diag(ones(length(S3),1)*Qfactor);
        STS = S3'*((Q.^2)\S3);
        C{idx} = STS;
        SAVE_deltaSOC(idx) =-trapz(t,DFN_out.Current*0.0744)/(3600*2.9);

        % Plot CURRENT & Save
% %         figure()
% %         plot(DFN_out.Current*0.0744)
% %         xlabel('Time')
% %         ylabel('Current')
% %         set(gca,'FontSize',16)
% %         save_name = [num2str(idx) '.png'];
% %         saveas(gcf,save_name)
% %         
% %         close all

        % Plot VOLTAGE & Save

% %         figure()
% %         plot(DFN_out.Volt)
% %         xlabel('Time')
% %         ylabel('Output Voltage')
% %         set(gca,'FontSize',16)
% %         save_name = [num2str(idx) '.png'];
% %         saveas(gcf,save_name)
% %         
% %         close all
        
        % Update index
        idx = idx +1;
        
    end
end


toc

%% BRIDGE
cost = transpose(cell2mat(cost));
%===========%===========%===========%===========%===========%===========

%% OED_CVX: D-optimality

% STS: Fisher Information Matrix (STS-form)
% cost: each experiment's cost (time)
% Profile: DFN output (struct)


num_inputs=length(C); % Total number of inputs to execute

m = 50; % The number of experiments to execute
Base_Budget_hr = 1*60*60; % Budget for an hour

MIN_m = 1;
MAX_m = 100; % Maximum number of experiments to execute
Nexp_exe = MIN_m : MAX_m;

obj_mat = 0;

MAX_HRs = 1*24;
Nexp_budget = 1:(10*24);

save_obj = zeros(length(Nexp_budget),length(Nexp_exe));

%% OED_CVX 3D SURFACE
% 
% for i=1:length(Nexp_budget)
%     fprintf('Iters: %d out of %d \n',i,length(Nexp_budget))
%     for j=1:length(Nexp_exe)
%         B = Base_Budget_hr * Nexp_budget(i);
%         m = Nexp_exe(j);
%         
%         cvx_begin quiet
%         cvx_solver SeDuMi 
% 
%         variable x(num_inputs);
% 
%         bigOne = ones(1,num_inputs);
%         expression M
% 
%         for k = 1 : num_inputs
%             M = M+x(k)*C{k};
%         end
% 
%         %maximize( log_det(M) ); % D-optimality [OK]
%         maximize(det_rootn(M)) % equivalent to D-optimality
%         % check details at http://ask.cvxr.com/t/log-determinant-base/104
%         
%         subject to
%             m*cost'*x <= B
%             bigOne*x == 1
%             0 <= x <= 1/m
% 
%         cvx_end
%         
%         if isnan(cvx_optval)
%             fprintf('OED_CVX failed at budget:%d [hrs], Nexp:%d \n',B,m);
%         end
%         
%         % SAVE & Recalculate objective
%         x_idx = round(m*x);
%         
%         for k= 1: num_inputs
%             obj_mat = obj_mat + x_idx(k)*C{k};
%         end
%         
%         save_obj(i,j) = det(obj_mat);
%         obj_mat =0;
%     end
% end
% 
% %
% figure()
% [X,Y] = meshgrid(1:length(Nexp_budget),1:length(Nexp_exe));
% surf(Nexp_budget(X),Nexp_exe(Y),save_obj');
% xlabel('Budget [Hrs]')
% ylabel('Number of Exp')
% zlabel('Objective (logdet)')
% set(gca,'FontSize',16)
% 
% % Choose optimal resource
% [Max_obj,I] = max(save_obj(:));
% [I_row, I_col] = ind2sub(size(save_obj),I); % I_row: budget (B), I_col: number of exp.(m)

%% DO OED-CVX-Q

B = Base_Budget_hr * 10;%Nexp_budget(I_row);
m = 5;%Nexp_exe(50);%Nexp_exe(I_col);

cvx_begin
    cvx_solver SeDuMi

    variable x(num_inputs);
    
    bigOne = ones(1,num_inputs);
    expression M
    
    for k = 1 : num_inputs
        M = M+x(k)*C{k};
    end
    
%     maximize( log_det(M) ); % D-optimality [OK]
%     minimize( log_det(inv(M)) ); % D-optimality [OK]
    maximize(det_rootn(M)) % equivalent to D-optimality
% check details at http://ask.cvxr.com/t/log-determinant-base/104
% https://github.com/cvxr/CVX/blob/master/functions/det_rootn.m
    
    subject to
        m*cost'*x <= B
        bigOne*x == 1
        0 <= x <= 1/m
        
cvx_end

Num_exp = round(m*x);
obj_mat = 0;
for k= 1: num_inputs
    obj_mat = obj_mat + Num_exp(k)*C{k};
end
fprintf('Determinant of STS: %e \n',det(obj_mat));

% if (logdet(obj_mat) ~= Max_obj)
%     fprintf('OED_CVX has error \n');
% end

% Plot Final OED_CVX output
figure()
bar(round(m*x))
xlabel('Experiment id #')
xlim([1,num_inputs])
ylabel('The Number of execution')
set(gca,'FontSize',16)
% title('OED CVX Q result')

% Plot Input library info.
for i=1:length(C)
    length_input(i) = length(Profile{i}.Current);
end

figure()
plot(length_input/3600,'LineWidth',2)
xlabel('Experiment identity number #')
xlim([1,num_inputs])
ylabel('Time length [hours]')
% title('TimeLength vs Profiles')
set(gca,'FontSize',16)


%% Result inspection

%%% Plot Input library info.
% for i=1:length(C)
%     det_input(i) = det(C{i});
% end
% figure()
% subplot(211)
% bar(det_input)
% xlabel('Exp. #')
% xlim([1,num_inputs])
% ylabel('Ind. Exp. Det')
% set(gca,'FontSize',16)
% subplot(212)
% bar(log(det_input))
% xlabel('Exp. #')
% xlim([1,num_inputs])
% ylabel('Ind. Exp. LogDet')
% set(gca,'FontSize',16)

%%% Plot Input library info.
% for i=1:length(C)
%     length_input(i) = length(Profile{i}.Current);
% end
% 
% figure()
% plot(length_input/3600)
% xlabel('Exp. #')
% xlim([1,num_inputs])
% ylabel('Individual Exp. length [hrs]')
% title('TimeLength vs Profiles')
% set(gca,'FontSize',16)

% % % Plot V0 in Input library
% % for i=1:length(Profile)
% %     V0_moeum(i) = Profile{i}.Volt(1);
% % end
% % 
% % figure()
% % plot(V0_moeum)
% % xlabel('Exp. #')
% % ylabel('Volt')
% % title('V0 vs Profiles')
%% OID
% run('OID_v1')

% figure;plot(currentprof)