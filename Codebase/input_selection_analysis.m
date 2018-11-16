%% Max Input Analysis
%%% By: Zach Gima 2018-11-1
%%% Code for analyzing inputs selected for sensitivity analysis

function input_selection_analysis(r,senspath)
    sens_inputs = struct([]);
    num_inputs = length(r);
    sens_index = cell(num_inputs,1);
    for i=1:num_inputs
        sens_index{i}=r(i).name;
        load(horzcat(senspath,sens_index{i}));

        %%% Need: t, I, V, V_0, T1, T2
        sens_inputs(i).t = t;
        sens_inputs(i).I = I;
        sens_inputs(i).T1 = alg_states.T1_sim - 273.15;
        sens_inputs(i).T2 = alg_states.T2_sim - 273.15;
        sens_inputs(i).T_a = T_a(1);
        sens_inputs(i).V_0 = V_0;
        sens_inputs(i).V = V;
        sens_inputs(i).maxT1 = max(sens_inputs(i).T1);
        sens_inputs(i).maxT2 = max(sens_inputs(i).T2);

    end    

    %% Input Stats
    V_0_count = zeros(1,4);
    T_a_count = zeros(1,4);
    exp_length_count = zeros(1,5); % [600 1350 2100 2850 3600]
    exp_type_count = zeros(1,3); % [charge discharge both]

    SOC_0_vec = [0.2 0.4 0.6 0.8];
    V_0_vec = [3.4562 3.5982 3.7689 3.9481]; %from ocv_soc_map.m
    T_a_vec = [0 20 40 60]; % just added 60 in even though it's not a temperature used
    exp_length_vec = [600 1350 2100 2850 3600];
    exp_type_cell = {'Charge', 'Discharge', 'Both'};

    for i = 1:num_inputs
        for jj = 1:4 % # of possible different T_amb and initial V
            % Check initial Voltage, T_amb, experiment length, and
            % charge/discharge
            V_0_count(jj) = V_0_count(jj) + (sens_inputs(i).V_0 == V_0_vec(jj));
            T_a_count(jj) = T_a_count(jj) + (sens_inputs(i).T_a == T_a_vec(jj));
        end

        %check input length
        if (length(sens_inputs(i).t) - 1 == exp_length_vec(1))
            exp_length_count(1) = exp_length_count(1) + 1;
        elseif (length(sens_inputs(i).t) - 1 == exp_length_vec(2))
            exp_length_count(2) = exp_length_count(2) + 1;
        elseif (length(sens_inputs(i).t) - 1 == exp_length_vec(3))
            exp_length_count(3) = exp_length_count(3) + 1;
        elseif (length(sens_inputs(i).t) - 1 == exp_length_vec(4))
            exp_length_count(4) = exp_length_count(4) + 1;
        else
            exp_length_count(5) = exp_length_count(5) + 1;
        end

        % check whether experiment charge/discharge/both
        if sum(sens_inputs(i).I) > 3600 % charge
            exp_type_count(1) = exp_type_count(1) + 1;
        elseif sum(sens_inputs(i).I) < -3600 % discharge
            exp_type_count(2) = exp_type_count(2) + 1;
        else %both case
            exp_type_count(3) = exp_type_count(3) + 1;
        end
    end



    %% Plot Everything
    fs = 25;

    folder_path  = 'Plots/';

    %%% Plot Selected Input Stats
    % Plot SOC_0 stats
    f = figure('Position', [100 100 900 700],'visible','off');
    bar(SOC_0_vec,V_0_count);
    xticks([0.2 0.4 0.6 0.8])
    yticks(0:1:max(V_0_count));
    title('Initial SOC for Max Sens. Inputs')
    xlabel('Initial SOC')
    ylabel('# Inputs Selected')
    set(gca,'Fontsize',fs)

    filename = strcat(folder_path,'SOC0_stats');
    saveas(f,filename,'fig')
    saveas(f,filename,'png')

    % Plot Temp stats
    f = figure('Position', [100 100 900 700],'visible','off');
    bar(T_a_vec(1:3),T_a_count(1:3));
    xticks([0 20 40])
    title('T_{amb} for Max Sens. Inputs')
    xlabel('T_{amb}')
    ylabel('# Inputs Selected')
    set(gca,'Fontsize',fs)

    filename = strcat(folder_path,'temp_stats');
    saveas(f,filename,'fig')
    saveas(f,filename,'png')

    % Plot exp length stats
    f = figure('Position', [100 100 900 700],'visible','off');
    bar(exp_length_vec,exp_length_count);
    xticks([600 1350 2100 2850 3600])
    title('Exp. Length for Max Sens. Inputs')
    xlabel('Experiment Length (s)')
    ylabel('# Inputs Selected')
    set(gca,'Fontsize',fs)

    filename = strcat(folder_path,'length_stats');
    saveas(f,filename,'fig')
    saveas(f,filename,'png')

    % Plot exp type stats
    f = figure('Position', [100 100 900 700],'visible','off');
    bar([1 2 3],exp_type_count);
    xticks([1 2 3])
    xticklabels({'Charge','Discharge','Both'})
    title('Exp. Type for Max Sens. Inputs')
    xlabel('Experiment Type')
    ylabel('# Inputs Selected')
    set(gca,'Fontsize',fs,'YTick',0:max(exp_type_count))

    filename = strcat(folder_path,'type_stats');
    saveas(f,filename,'fig')
    saveas(f,filename,'png')


    %%% Plot all of the individual Inputs with voltage and temp. responses
    for pp = 1:num_inputs

        input_num_str = sens_index{pp}(1:end-4);

        % Plot Input (t vs. I)
        f_I = figure('Position', [100 100 900 700],'visible','off');
        plot(sens_inputs(pp).t,sens_inputs(pp).I,'r','LineWidth',3)

        title_text = strcat('Current, Input Profile: #',input_num_str);
        title(title_text)
        xlabel('Time (s)')
        ylabel('Input Current (A)')
        set(gca,'Fontsize',fs) 

        filename_I = strcat(folder_path,input_num_str,'_I');
        saveas(f_I,filename_I,'fig')
        saveas(f_I,filename_I,'png') 

        f_V = figure('Position', [100 100 900 700],'visible','off');
        plot(sens_inputs(pp).t,sens_inputs(pp).V,'g','LineWidth',3)

        title_text = strcat('Voltage Response, Input Profile: #',input_num_str);
        title(title_text)
        xlabel('Time (s)')
        ylabel('Voltage (V)')
        set(gca,'Fontsize',fs)

        filename_V = strcat(folder_path,input_num_str,'_V');
        saveas(f_V,filename_V,'fig')
        saveas(f_V,filename_V,'png')  

        % Plot Temperature evolution (t vs. T)
        f_T = figure('Position', [100 100 900 700],'visible','off');
        hold on
        plot(sens_inputs(pp).t,sens_inputs(pp).T1,'--b','LineWidth',3)
        plot(sens_inputs(pp).t,sens_inputs(pp).T2,'-r','LineWidth',3)
        hold off

        legend('T1','T2')
        title_text = strcat('Temperature, Input Profile: #',input_num_str);
        title(title_text)
        xlabel('Time (s)')
        ylabel('Temperature (C)')
        set(gca,'Fontsize',fs) 

        filename_T = strcat(folder_path,input_num_str,'_T');
        saveas(f_T,filename_T,'fig')
        saveas(f_T,filename_T,'png') 
    end
end