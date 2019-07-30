clear all
clc
fs = 25;
io_folder = '/Users/ztakeo/Documents/GitHub/OED_ParamID/ID_results/';

% load experimental + other data
load('InputLibrary/ValidationCycles/validation_cycles_fmincon.mat');

% load simulated val. cycles for nominal parameter set
load('validation_results_nominal.mat','V_LM_CELL_sim');
V_LM_CELL_nominal = V_LM_CELL_sim(17:end)';
clear V_LM_CELL_sim

% load simulated val. cycles for ID'ed parameter set
load('fmincon_params.mat');

num_exp = length(V_LM_CELL);

for ii = 1:num_exp
    profile_name = exp_num{ii}(1:end-4);

    % Plot voltage fit
    figure('Position', [100 100 900 700])
    if ii < 7
        plot(Time_exp{ii},V_LM_CELL{ii},Time_exp{ii},V_LM_CELL_sim{ii},Time_exp{ii},V_LM_CELL_nominal{ii},'LineWidth',2)
        % Calculate RMSE
        rmse_save_nominal(ii,1) = rmse(V_LM_CELL{ii},V_LM_CELL_nominal{ii});
        rmse_save_ID(ii,1) = rmse(V_LM_CELL{ii},V_LM_CELL_sim{ii});
    else
        idx_end = length(V_LM_CELL_nominal{ii});
        plot(Time_exp{ii}(1:idx_end),V_LM_CELL{ii}(1:idx_end),Time_exp{ii}(1:idx_end),V_LM_CELL_sim{ii}(1:idx_end),Time_exp{ii}(1:idx_end),V_LM_CELL_nominal{ii},'LineWidth',2)
        % Calculate RMSE
        rmse_save_nominal(ii,1) = rmse(V_LM_CELL{ii}(1:idx_end),V_LM_CELL_nominal{ii});
        rmse_save_ID(ii,1) = rmse(V_LM_CELL{ii},V_LM_CELL_sim{ii});
    end
    %     legend('DFN','Identified Params')
    legend('Experimental Data','Identified Params', 'Nominal Params','Location','northwest','NumColumns',3)
    xlabel('Time (s)','Fontsize',fs)
    ylabel('Voltage (V)','Fontsize',fs)
%     title('Voltage Fit')
    set(gca,'FontSize',fs)
    box on
    grid on

%     % Calculate RMSE
%     rmse_save_nominal(ii,1) = rmse(V_LM_CELL{ii},V_LM_CELL_nominal{ii});
%     rmse_save_ID(ii,1) = rmse(V_LM_CELL{ii},V_LM_CELL_sim{ii});
    
    % Save figures
    savefig(strcat(io_folder,profile_name,'_voltage_fit.fig'))
    print(strcat(io_folder,profile_name,'_voltage_fit'),'-dpng')
    
    close all
end