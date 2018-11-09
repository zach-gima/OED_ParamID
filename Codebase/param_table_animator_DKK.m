%% Function for plotting the true parameters and the estimated parameters with confidence intervals
% Adapted from param_table_figure.m, created by Saehong Park
% By: Zach Gima 2018-3-20

function param_table_animator_DKK(ci95,theta_0,theta_0_true,truth_param,bounds)
    %% Fig setting
    fs = 25;

    %% parameter numbering
    % Dynamic parameter only numbering
    %1. D_s_n
    %2. D_s_p
    %3. R_s_n
    %4. R_s_p
    %5. sig_n
    %6. sig_p
    %7. D_e
    %8. eps_e_n
    %9. epe_e_s
    %10. eps_e_p
    %11. Kappa
    %12. t_plus
    %13. d_activity
    %14. k_n0
    %15. k_p0
    %16. R_f_n
    %17. R_f_p
    %18. ce0

    % params = {'$D_s^{^{\_}}$','$D_s^+$','$R_s^-$','$R_s^+$',...
    %     '$\sigma^{^{\_}}$','$\sigma^{+}$','$D_e(\cdot)$','$\varepsilon_e^{^{\_}}$','$\varepsilon_e^{sep}$','$\varepsilon_e^+$',...
    %     '$\kappa(\cdot)$','$t_{c}^{0}$','$\frac{d \ln f_{c/a}}{d \ln c_e}(\cdot)$','$k^{^{\_}}$','$k^+$','$R_f^{^{\_}}$','$R_f^+$','$c_{e_0}$'};

    %%%group4->group3->group2->group1
    % params = {'$t_{c}^{0}$','$\varepsilon_e^{sep}$','$\sigma^{+}$','$\sigma^{^{\_}}$','$k^+$','$k^{^{\_}}$','$c_{e_0}$',...
    %             '$\varepsilon_e^+$','$R_f^+$','$R_f^{^{\_}}$','$\frac{d \ln f_{c/a}}{d \ln c_e}(\cdot)$','$D_e(\cdot)$',...
    %             '$\kappa(\cdot)$','$\varepsilon_e^{^{\_}}$','$D_s^+$','$D_s^{^{\_}}$','$R_s^+$','$R_s^-$'};
    % params_idx_from_original = [12;9;6;5;15;14;18;10;17;16;13;7;11;8;2;1;4;3]; % params index from params original        

    %%%group4->group3->group2->group1 align with paper frame.
    params = {'$t_{c}^{0}$','$\sigma^{^{\_}}$','$\sigma^{+}$','$k^+$','$\varepsilon_e^{sep}$','$k^{^{\_}}$','$R_f^+$','$c_{e_0}$',...
                '$\varepsilon_e^+$','$R_f^{^{\_}}$','$\frac{d \ln f_{c/a}}{d \ln c_e}(\cdot)$','$D_s^{^{\_}}$','$D_e(\cdot)$',...
                '$\varepsilon_e^{^{\_}}$','$\kappa(\cdot)$','$D_s^+$','$R_s^-$','$R_s^+$'};
    params_idx_from_original = [12;5;6;15;9;14;17;18;10;16;13;1;7;8;11;2;3;4];        

    selection_vector = [1;1;1;1;0;0;1;1;1;1;1;1;1;1;1;1;1;1;1;0;1]; %G4
    sel_k = find(selection_vector(:,end));


    %% parameter values
    % Load truth parameter values, truth_param
    Selected_truthParam = truth_param(sel_k);             
    
    % Load estimated parameters
    Final_param = theta_0;
    Selected_params = Final_param(sel_k);
    
    % Load Initial Parameter Guess [ZTG Change]
%     run param/params_nominal
    Selected_initial_params = theta_0_true(sel_k);                   

    %Load confidence intervals
    %Final_CI = ci95;

    %Final_CI = Final_CI + Final_param;        
    %Selected_CI = Final_CI(sel_k);  
    
    %% Parameter bound  
    normalized_truth_bar = origin_to_norm(Selected_truthParam,bounds,selection_vector);
    normalized_theta_bar = origin_to_norm(Selected_params,bounds,selection_vector);                     
    %normalized_CI_bar = origin_to_norm(Selected_CI,bounds,selection_vector);

%     % NOTE: IF USING origin_to_norm_MinMax(...), WILL HAVE TO REWRITE IT
%     SO THAT IT ACCOMODATES THE NEW PARAMS_BOUNDS NOTATION
%     MinMaxNormalized_truth_bar = origin_to_norm_MinMax(Selected_truthParam,bounds,selection_vector);
%     MinMaxNormalized_theta_bar = origin_to_norm_MinMax(Selected_params,bounds,selection_vector);                     
%     MinMaxNormalized_CI_bar = origin_to_norm_MinMax(Selected_CI,bounds,selection_vector);
    
    %[ZTG Change]
    normalized_initial_theta_bar = origin_to_norm(Selected_initial_params,bounds,selection_vector);                     
%     MinMaxNormalized_initial_theta_bar = origin_to_norm_MinMax(Selected_initial_params,bounds,selection_vector);  

    %% Plot %1

    % convert index

    x_value = normalized_theta_bar(params_idx_from_original)';
    y_value = 1:18;
    truth_value = normalized_truth_bar(params_idx_from_original)';
    %err = normalized_CI_bar(params_idx_from_original)';
    %err = x_value - err;
%     err(3) = -2; % This is because different normalization scheme. i.e, log-normal. So manually change it.
    %err(2) = -2; % This is because different normalization scheme. i.e, log-normal. So manually change it.
    Np = 18;
    
    %[ZTG Change]
    initial_value = normalized_initial_theta_bar(params_idx_from_original);

    x_origin = 0:0.1:1;
    y_dotted = [1*ones(1,length(x_origin));
                2*ones(1,length(x_origin));
                3*ones(1,length(x_origin));
                4*ones(1,length(x_origin));
                5*ones(1,length(x_origin));
                6*ones(1,length(x_origin));
                7*ones(1,length(x_origin));
                8*ones(1,length(x_origin));
                9*ones(1,length(x_origin));
                10*ones(1,length(x_origin));
                11*ones(1,length(x_origin));
                12*ones(1,length(x_origin));
                13*ones(1,length(x_origin));
                14*ones(1,length(x_origin));
                15*ones(1,length(x_origin));
                16*ones(1,length(x_origin));
                17*ones(1,length(x_origin));
                18*ones(1,length(x_origin))];
    
    %if (runonce)
    figure(2); clf;
%     set(gcf,'Position',[234 3 564 695],'PaperPositionMode','auto');
    set(gcf,'Position',[100 100 900 700],'PaperPositionMode','auto');

    plot(x_origin,y_dotted,'k--') % BaseLine
    set(gca,'YTick',1:Np);
    set(gca,'Position',[0.2 0.1 0.75 0.85]);
    [hx,hy] = format_ticks(gca,' ', params(y_value),0:0.1:1,[],0,0,0.05,'FontSize',fs,'FontWeight','Bold');
    xlabel('$\bar{\theta}$','Interpreter','Latex','FontSize',16);
    hold on
    e=plot(x_value,y_value,'LineStyle','none');%,'horizontal','o');
    ylim([0 20])
    e.Marker = '*';
    e.MarkerSize = 10;
    e.Color = 'red';
    e.LineWidth = 2;
    %runonce = false;
    

    %end
    hold on;
    
       % Plot truth and initial values
    for ii = 1:length(truth_value)
        plot(truth_value(ii),ii,'bo','LineWidth',2,'MarkerSize',15) % Truth value.
        plot(initial_value(ii),ii,'k^','LineWidth',2,'MarkerSize',15) % Initial values.
    end

   %,'Confidence interval');
   h = zeros(3, 1);
    h(1) = plot(NaN,NaN,'ob','LineWidth',2);
    h(2) = errorbar(NaN,NaN,'r*','LineWidth',2);
    h(3) = plot(NaN,NaN,'k^','LineWidth',2);
    %h(4) = plot(NaN,NaN,'r-','LineWidth',2);
    legend(h,'Truth value','Final estimate','Initial estimate');
end