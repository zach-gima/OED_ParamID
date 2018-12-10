%% Function for plotting the true parameters and the estimated parameters with confidence intervals
% Adapted from param_table_figure.m, created by Saehong Park
% By: Zach Gima 2018-3-20

%%%%%%% Note param_table_plotter should always use final group's selection_vector
function [g] = param_table_plotter_ZTG(mode,selection_vector,ci95,Final_param,theta_0_true,truth_param,bounds)
    fs = 25;
    sel_k = find(selection_vector);

    Np = length(sel_k); % Number of ID'ed params (less than Np_full)
    Np_full = length(theta_0_true); % Number of possible params to identify

    %% parameter numbering
    % ZTG Updated: 2018-12-10
    % Blank means not selected to identify
    % Eq: equilibrium structure
    
    %   (G1) 1  : D_s_n       => p.D_s_n0
    %   (G1) 2  : D_s_p       => p.D_s_p0
    %   (G1) 3  : R_s_n       => p.R_s_n
    %   (G1) 4  : R_s_p       => p.R_s_p
    %   (Eq) 5  : epsilon_s_n => p.epsilon_s_n
    %   (Eq) 6  : epsilon_s_p => p.epsilon_s_p
    %   () 7  : sig_n       => p.sig_n
    %   () 8  : sig_p       => p.sig_p
    %   (G1) 9  : D_e         => p.ElecFactorD --> multiply factor
    %   () 10 : epsilon_e_n => p.epsilon_e_n
    %   () 11 : epsilon_e_s => p.epsilon_e_s
    %   () 12 : epsilon_e_p => p.epsilon_e_p
    %   () 13 : kappa       => p.ElecFactorK --> multiply factor
    %   () 14 : t_plus      => p.t_plus
    %   (G1) 15 : dactivity   => p.ElecFactorDA --> multiply factor
    %   () 16 : k_n0        => p.k_n0
    %   () 17 : k_p0        => p.k_p0
    %   () 18 : R_f_n       => p.R_f_n
    %   () 19 : R_f_p       => p.R_f_p
    %   (Eq) 20 : n_Li_s      => p.n_Li_s
    %   () 21 : c_e0        => p.c_e
    %   () 22 : E.Dsn        => p.E.Dsn
    %   () 23 : E.Dsp        => p.E.Dsp
    %   () 24 : E.kn        => p.E.kn
    %   () 25 : E.kp        => p.E.kp

    % Full 25 parameters that we could identify
    params_full = {'$D_s^{^{\_}}$','$D_s^+$','$R_s^-$','$R_s^+$','$\varepsilon_s^{^{\_}}$','$\varepsilon_s^+$',...
        '$\sigma^{^{\_}}$','$\sigma^{+}$','$D_e(\cdot)$','$\varepsilon_e^{^{\_}}$','$\varepsilon_e^{sep}$','$\varepsilon_e^+$',...
        '$\kappa(\cdot)$','$t_{c}^{0}$','$\frac{d \ln f_{c/a}}{d \ln c_e}(\cdot)$','$k^{^{\_}}$','$k^+$','$R_f^{^{\_}}$',...
        '$R_f^+$','$n_{Li_s}$','$c_{e_0}$','$E.D_s^-$','$E.D_s^+$','$E.k_n$','$E.k_p$'};
   
    % Parameters remaining to identify after removing eq. struct. params as well as params eliminated during
    % clustering and sens. threshold analysis
    %%%Ordered from group2->group1 in terms of sens. magnitude i.e. last
    %%%param in cell is highest sens.
    params_ID = {'$E.k_n$','$E.D_s^-$','$\varepsilon_e^{sep}$','$c_{e_0}$','$\varepsilon_e^+$','$E.D_s^+$','$R_f^{^{\_}}$',...
            '$\frac{d \ln f_{c/a}}{d \ln c_e}(\cdot)$','$D_s^{^{\_}}$','$D_e(\cdot)$','$D_s^+$','$R_s^-$','$R_s^+$'};
    
    params_idx_from_original = zeros(Np,1);
    
    for ii = 1:Np
        params_idx_from_original(ii) = find(strcmp(params_ID{ii},params_full) == 1);
    end

    %% parameter values

    % Load truth parameter values, truth_param
    Selected_truthParam = truth_param(sel_k);             
    
    % Load estimated parameters
    Selected_params = Final_param(sel_k);
    
    % Load Initial Parameter Guess [ZTG Change]
%     run param/params_nominal
    Selected_initial_params = theta_0_true(sel_k);                  

    % Load confidence intervals
    Final_CI = ci95;

    Final_CI = Final_CI + Final_param;        
    Selected_CI = Final_CI(sel_k);  
    
    %% Parameter bound
    normalized_truth_bar = zeros(Np_full,1);
    normalized_theta_bar = zeros(Np_full,1);
    normalized_CI_bar = zeros(Np_full,1);
    normalized_initial_theta_bar = zeros(Np_full,1);
    
    % These variables are Num. Param x 1 length;
    normalized_truth_bar(sel_k) = origin_to_norm('param',Selected_truthParam,bounds,selection_vector);
    normalized_theta_bar(sel_k) = origin_to_norm('param',Selected_params,bounds,selection_vector);                     
    normalized_CI_bar(sel_k) = origin_to_norm('sens',Selected_CI,bounds,selection_vector);

%     % NOTE: IF USING origin_to_norm_MinMax(...), WILL HAVE TO REWRITE IT
%     SO THAT IT ACCOMODATES THE NEW PARAMS_BOUNDS NOTATION
%     MinMaxNormalized_truth_bar = origin_to_norm_MinMax(Selected_truthParam,bounds,selection_vector);
%     MinMaxNormalized_theta_bar = origin_to_norm_MinMax(Selected_params,bounds,selection_vector);                     
%     MinMaxNormalized_CI_bar = origin_to_norm_MinMax(Selected_CI,bounds,selection_vector);
    
    %[ZTG Change]
    normalized_initial_theta_bar(sel_k) = origin_to_norm('param',Selected_initial_params,bounds,selection_vector);                     
%     MinMaxNormalized_initial_theta_bar = origin_to_norm_MinMax(Selected_initial_params,bounds,selection_vector);  

    %% Plot %1

    % convert index
    truth_value = normalized_truth_bar(params_idx_from_original)';
    x_value = normalized_theta_bar(params_idx_from_original)';
    err = normalized_CI_bar(params_idx_from_original)';
    err = x_value - err;
%     err(3) = -2; % This is because different normalization scheme. i.e, log-normal. So manually change it.
%     err(2) = -2; % This is because different normalization scheme. i.e, log-normal. So manually change it.
    initial_value = normalized_initial_theta_bar(params_idx_from_original);

    % Dartboard plotting variables
    x_origin = 0:0.1:1;
    y_value = 1:Np;
    y_dotted = zeros(Np,length(x_origin));

    for ii = 1:Np
        y_dotted(ii,:) = ii*ones(1,length(x_origin));
    end
    
    if (strcmp(mode,'static') == 1)
        figure(); clf;
    elseif (strcmp(mode,'animate') == 1)
        g = figure(2);clf;
    else
        error('Specify either static plot (static) or GIF (animate)');
    end
    
%     set(gcf,'Position',[234 3 564 695],'PaperPositionMode','auto');
    set(gcf,'Position',[100 100 900 700],'PaperPositionMode','auto');

    plot(x_origin,y_dotted,'k--') % BaseLine
    set(gca,'YTick',1:Np);
    set(gca,'Position',[0.2 0.1 0.75 0.85]);
    [~,~] = format_ticks(gca,' ', params_ID(y_value),0:0.1:1,[],0,0,0.05,'FontSize',fs,'FontWeight','Bold');
    xlabel('$\bar{\theta}$','Interpreter','Latex','FontSize',16);
    hold on
    ylim([0 Np])
    
    if (strcmp(mode,'static') == 1)
        e=errorbar(x_value,y_value,err,'horizontal','o');
        e.Marker = '*';
        e.MarkerSize = 15;
        e.Color = 'red';
        e.CapSize = 15;
        e.LineWidth = 2;
    elseif (strcmp(mode,'animate') == 1)
        e=plot(x_value,y_value,'LineStyle','none');%,'horizontal','o');
        e.Marker = '*';
        e.MarkerSize = 10;
        e.Color = 'red';
        e.LineWidth = 2;
    end
    
    hold on;
    
    % Plot truth and initial values
    for ii = 1:length(truth_value)
        plot(truth_value(ii),ii,'bo','LineWidth',2,'MarkerSize',15) % Truth value.
        plot(initial_value(ii),ii,'k^','LineWidth',2,'MarkerSize',15) % Initial values.
    end
    
    h(1) = plot(NaN,NaN,'ob','LineWidth',2);
    h(2) = errorbar(NaN,NaN,'r*','LineWidth',2);
    h(3) = plot(NaN,NaN,'k^','LineWidth',2);
    if (strcmp(mode,'static') == 1)
        h(4) = plot(NaN,NaN,'r-','LineWidth',2);
        legend(h,{'Truth value','Final estimate','Initial estimate','Confidence interval'},'Fontsize',15);
    elseif (strcmp(mode,'animate') == 1)
        legend(h,{'Truth value','Final estimate','Initial estimate'},'Fontsize',15);

    end
%     set(gca,'FontSize',fs)

    hold off;
    
end