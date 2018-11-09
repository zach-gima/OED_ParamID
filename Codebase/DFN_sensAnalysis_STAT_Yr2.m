%% Plot Sensitivity Analysis Results of DFN
%   Created February 23, 2014 by Scott Moura
%   Modified May 11, 2017 by Saehong Park

clear;
close all;
clc;

fs = 18;

% addpath DFN

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
%   19 : E.Dsn       => p.E.Dsn
%   20 : E.Dsp       => p.E.Dsp
%   21 : E.kn        => p.E.kn
%   22 : E.kp        => p.E.kp

params = {'$D_s^{^{\_}}$','$D_s^+$','$R_s^-$','$R_s^+$',...
    '$\sigma^{^{\_}}$','$\sigma^{+}$','$D_e(\cdot)$','$\varepsilon_e^{^{\_}}$','$\varepsilon_e^{sep}$','$\varepsilon_e^+$',...
    '$\kappa(\cdot)$','$t_{c}^{0}$','$\frac{d \ln f_{c/a}}{d \ln c_e}(\cdot)$','$k^{^{\_}}$','$k^+$','$R_f^{^{\_}}$','$R_f^+$','$c_{e_0}$',...
    '$E.D_s^-$','$E.D_s^+$','$E.k_n$','$E.k_p$'};

run param/params_NCA

%% Normalization

% % % % Log normalization
% % % % Multiply (1/log10(exp(1))) * log(\theta_max / \theta_min) * theta_0

% % % % minmax normalization
% % % % Multiply (\theta_max - \theta_min)

% % % % nominal range normalization
% % % % Multiply norminal param (\theta_0)

% % % % Parameter Min/Max values

Dsn_max = 1.05E-12; % LCO, since NCA literature is too small
Dsn_min = 2.25E-16; % LCO, ""
Dsp_max = 1.00E-12; % LCO, ""
Dsp_min = 2.00E-16; % LCO, ""

Rsn_min = 1.0e-06;
Rsn_max = 100e-06;
Rsp_min = 1.0e-06;
Rsp_max = 100e-06;

Sig_n_max = 500; % NCA, times 5 by nominal param
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

EDsn_min = 22e3;
EDsn_max = 48.9e3;
EDsp_min = 20e3;
EDsp_max = 80.6e3;

Ekn_min = 37.48e3;
Ekn_max = 67.5e3;
Ekp_min = 30e3;
Ekp_max = 42.6e3;

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
    EDsn_max - EDsn_min; ...
    EDsp_max - EDsp_min; ...
    Ekn_max - Ekn_min; ...
    Ekp_max - Ekp_min];
normalize_param = normalize_param' ; % row maxtrix

%% Load Data;
% Load SensG*_*.mat

Save_S = {};
Save_S_GS = {};
Save_Snorm = {};
SS_idx = 1;
%%%%%%%%%%%%%%%%%%  %%%%%%%%%%%%%%%%%%%%%%
% Load Sensitivities

% fn= 'All\';
fn= 'SensResults/All/';
disp(['Loaded Sensitivity directory:  ' fn]);

for i = 1:1615
    
    data_path = [fn num2str(i) '.mat'];
    try
        load(data_path); % stuct file: sens
    catch
        continue;
        %         fprintf(['Missing: ' data_path '\n']);
        %         sens.S3 = 1;
        %         sens.DFNout = 1;
        %         sens.COST = 1e6;
        %         sens.Err_flag = 1;
    end
    S3 = S_CASADI;
    
    for k = 1:length(S_CASADI)
        S3(k,:) = normalize_param .* S3(k,:);
    end
    
    %     S3 = normalize_param .* S3;
    
    %     DFN_out = sens.DFNout;
    %     Cost = sens.COST;
    %     t= linspace(1,Cost,Cost);
    Np = size(S3,2); % Number of parameters
    STS = S3'*S3;
    %     Err_flag = sens.Err_flag;
    
    if max(alg_states.T2_sim) >= 273.15+45
        fprintf('Input %d violates temperature condition, max tempearture is %f (T2) \n', i, max(alg_states.T2_sim - 273));
        continue;
    end
    
    Snorm.volt.unsort = zeros(Np,1);
    for idx = 1: Np
        Snorm.volt.unsort(idx) = norm(S3(:,idx)); % Non-orthnormalized version, not used here, but as comparison
    end
    
    % GranSchmidt Orthonormalization (GS)
    [Q,R,E] = qr(S3,0);
    D = diag(R);
    S_rank_log_ortho = flipud(log10(abs(D)));
    E_original = fliplr(E); % order of GS
    for j = 1:Np
        Unsort_S_rank_log_ortho(j) = S_rank_log_ortho(find(E_original == j));
    end
    
    % Save sensitivity
    Save_S{SS_idx} = Snorm.volt.unsort;
    Save_S_GS{SS_idx} = Unsort_S_rank_log_ortho;
    SS_idx = SS_idx+1;
    
    [Snorm.volt.sort, Snorm.volt.ind] = sort(Snorm.volt.unsort,1,'descend');
    S_rank_log = log10(Snorm.volt.sort(end:-1:1));
    S_rank_log_unsort = log10(Snorm.volt.unsort);
    
end


%% Sensitivity statistics for GramSchmidt(GS) version

disp('Plot sensitivity statistics?')

% Get random cell array with n cells. Here n = 150
n = size(Save_S_GS,2);

% Get size of matrices in cell.
matSize = size(Save_S_GS{1},1);
B = reshape(cell2mat(Save_S_GS),matSize,[],n);

% Sum 3D matrix along 3rd dimension
C = sum(B,3)./n;

% Sorting and Plot.
% Be careful to the order.
[GS_Snorm_sort, GS_Snorm_ind] = sort(C,2,'descend');
GS_rank_log = GS_Snorm_sort(end:-1:1);

figure(1); clf;
set(gcf,'Position',[234     3   564   695],'PaperPositionMode','auto');

barh(GS_rank_log);
set(gca,'YTick',1:Np);
set(gca,'Position',[0.2 0.1 0.75 0.85])
%set(gca, 'YTickLabel', params);
[hx,hy] = format_ticks(gca,' ',params(GS_Snorm_ind(end:-1:1)),-15:1:0,[],0,0,0.05,'FontSize',fs,'FontWeight','Bold');
set(gca,'FontSize',fs);
xlabel('Sensitivity Magnitude [log scale]','FontSize',fs)
% title('\bf Orthonormalized Sensitivity of Voltage','FontSize',fs+2);









