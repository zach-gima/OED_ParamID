%% plot stuff

% Note: this function written to work in Matlab 2016b and
% newer; issue w/ errorbar class

function Param_ID_plot(truth_param,theta_0_true,sel_k,paramID_out,ci95_full,t_paramID,rmse_final,output_folder,LM_options,bounds,alg_states,selection_vector)

    % Parse L-M Conditions
    param_exit_thresh = LM_options.exit_cond(1);
    chi_sq_exit_thresh = LM_options.exit_cond(2);
        
    fs = 25;
    iter = 1:length(paramID_out.save_RMSE);

    %% ParamID Dartboard
    park0 = paramID_out.save_param_org(:,end); %loads the identified parameters. 

    Final_param = zeros(length(theta_0_true),1);
    Final_param(sel_k) = park0;
    
    %%% Static Figure
    param_table_plotter_ZTG('static',selection_vector(:,end),ci95_full,Final_param,theta_0_true,truth_param,bounds);

    savefig(strcat(output_folder,'param_estimates.fig'))
    print(strcat(output_folder,'param_estimates'),'-dpng')

    %% GIF
%     animate_DKK(output_folder,selection_vector,truth_param)

    %% Plot Individual Parameter Evolution
    
    %% Plot LM logic
    if length(paramID_out.LM_logic) ~= iter(end)
        paramID_out.LM_logic = vertcat(paramID_out.LM_logic,0);
    end
    
    figure('Position', [100 100 900 700])
    plot(iter, paramID_out.LM_logic,'-o','LineWidth',3,'MarkerSize',10,'MarkerEdgeColor','k');
%     title('RMSE vs. Wall Clock time')
    xlabel('Iteration')
    title('Algorithm Action')
    ylim([-1 1])
    xlim([1,iter(end)]);
    
    actions = {'Decrease \lambda'; 'Recalculate Sens'; 'Increase \lambda'};
    set(gca,'Fontsize',fs,'ytick',[-1:1],'yticklabel',actions,'XTick',[1 : 1 : iter(end)])
    box on
    grid on
    savefig(strcat(output_folder,'LM_logic.fig'))
    print(strcat(output_folder,'LM_logic'),'-dpng')
    
    %% Plot Cost Function
    figure('Position', [100 100 900 700])
    plot(iter, paramID_out.save_chi_sq,'-o','LineWidth',3,'MarkerSize',10,'MarkerEdgeColor','k');
    title('Cost Function vs. Iteration')
    xlabel('Iteration')
    ylabel('Cost Function')
    %             ylim([0 0.025])
    xlim([1,iter(end)]);
    set(gca,'Fontsize',fs,'XTick',[1 : 1 : iter(end)])
    box on
    grid on

    savefig(strcat(output_folder,'cost_fcn.fig'))
    print(strcat(output_folder,'cost_fcn'),'-dpng')
    
%     %% Plot algebraic states rmse
%     %%% need to fix this by pulling only some of the values from the cells
%     %%% (dimension mismatch currently)
%     Num_inputs = length(paramID_out.Time_exp);
%     cssn_sim = cell(Num_inputs,1);
%     cssp_sim = cell(Num_inputs,1);
%     etan_sim = cell(Num_inputs,1);
%     etap_sim = cell(Num_inputs,1);
%     
%     for idx = 1:Num_inputs
%        % Parse specific algebraic states of interest
%        cssn_sim{idx} = alg_states{idx}.cssn_sim;
%        cssp_sim{idx} = alg_states{idx}.cssp_sim;
%        etan_sim{idx} = alg_states{idx}.etan_sim;
%        etap_sim{idx} = alg_states{idx}.etap_sim;
%     end
    
    %% plot RMSE evolution vs. time
    % use vertcat for adding in very initial RMSE
    figure('Position', [100 100 900 700])
    plot(t_paramID/3600, rmse_final,'-o','LineWidth',3,'MarkerSize',10,'MarkerEdgeColor','k');
    title('RMSE vs. Wall Clock time')
    xlabel('Wall Clock time (Hr)')
    ylabel('RMSE (V)')
    set(gca,'Fontsize',fs)
    box on
    grid on

    savefig(strcat(output_folder,'rmse_evolution.fig'))
    print(strcat(output_folder,'rmse_evolution'),'-dpng')

    %% Plot Normalized Distance between True and Estimated Parameters
    norm_param_dist = zeros(size(paramID_out.save_param_nmz,2),1);
    norm_truth_param = origin_to_norm('param',truth_param(sel_k),bounds,selection_vector(:,end));
    
    for ii = 1:size(paramID_out.save_param_nmz,2)
        norm_param_dist(ii) = norm(norm_truth_param-paramID_out.save_param_nmz(:,ii),2);
    end
    
    figure('Position', [100 100 900 700])
    plot(iter, norm_param_dist,'-o','LineWidth',3,'MarkerSize',10,'MarkerEdgeColor','k');
    title('Normalized Parameter Guess Error')
    xlabel('Iteration','FontSize',fs)
    ylabel('$|| \theta^* - \hat\theta ||$','Interpreter','Latex','FontSize',fs)  
    xlim([1,iter(end)]);
    set(gca,'Fontsize',fs,'XTick',[1 : 1 : iter(end)])

    box on
    grid on

    savefig(strcat(output_folder,'norm_param_dist.fig'))
    print(strcat(output_folder,'norm_param_dist'),'-dpng')
    
    %% plot CDF
    figure('Position', [100 100 900 700])
    h(1) = cdfplot(abs(paramID_out.save_y_minus_yfit(:,end)));
    set(gca,'Fontsize',fs)
    set(h,'LineWidth',3)
    xlabel('Voltage Error (V)')
    
    savefig(strcat(output_folder,'cdf_plot.fig'))
    print(strcat(output_folder,'cdf_plot'),'-dpng')
    
    %% Voltage Data Fit

    %%% Plot initial voltage fit (initial parameter guess vs. t and truth
    %%% voltage vs. t)
    
    
    %%% Y_dat & Y_sim
    figure('Position', [100 100 900 700])
    t = 1:length(paramID_out.y_dat);
    plot(t,paramID_out.y_dat,t,paramID_out.save_y_sim(:,end),'LineWidth',2)
%     legend('DFN','Identified Params')
    legend('Experimental Data','Identified Params')
    xlabel('Time (s)','Fontsize',fs)
    ylabel('Voltage (V)','Fontsize',fs)
    title('Voltage Fit')
    set(gca,'FontSize',fs)
    box on
    grid on

    savefig(strcat(output_folder,'voltage_fit.fig'))
    print(strcat(output_folder,'voltage_fit'),'-dpng')
    
    %%% Voltage RMSE
    figure('Position', [100 100 900 700])
    plot(iter,paramID_out.save_RMSE,'-o','LineWidth',3, 'MarkerSize',10,'MarkerEdgeColor','k')
    xlabel('Iteration','Fontsize',fs)
    ylabel('Voltage RMSE (V)','Fontsize',fs)
    title('Voltage RMSE vs. Iter')

    xlim([1,iter(end)]);
    set(gca,'Fontsize',fs,'XTick',[1 : 1 : iter(end)])    
    box on
    grid on
    
    savefig(strcat(output_folder,'voltage_rmse.fig'))
    print(strcat(output_folder,'voltage_rmse'),'-dpng')
    
    %% Exit Condition Evolution

    figure('Position', [100 100 900 700])
    hold on
    plot(iter,paramID_out.save_param_exit,'-o','LineWidth',3,'MarkerSize',10,'MarkerEdgeColor','k')
    plot(iter,param_exit_thresh*ones(size(paramID_out.save_param_exit)),'k--','LineWidth',3) %plot parameter convergence exit criterion

    plot(iter,paramID_out.save_chi_sq_exit,'-+','LineWidth',3,'MarkerSize',10,'MarkerEdgeColor','k')
    plot(iter,chi_sq_exit_thresh*ones(size(paramID_out.save_chi_sq_exit)),'k:','LineWidth',3) %plot parameter convergence exit criterion

    xlabel('Iteration','FontSize',fs)
    ylabel('Normalized Exit Conditions','FontSize',fs)
    title('Exit Condition Evolutions, G4','FontSize',fs)
    legend('Parameter Convergence', 'Parameter Exit Threshold','Cost Function Convergence','Cost Function Exit Threshold')
    
    xlim([1,iter(end)]);
    set(gca,'FontSize',fs,'yscale','log','XTick',[1 : 1 : iter(end)])
    box on
    grid on
    hold off
    
    savefig(strcat(output_folder,'exit_conditions.fig'))
    print(strcat(output_folder,'exit_conditions'),'-dpng')
end