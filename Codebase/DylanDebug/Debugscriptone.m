   clear e_idx rand_idx25
   goupsize = 2
for i = 1:5
    W = [1;1;1;1;0;0;0;0;1;0;1;1;0;0;1;0;0;1;0;0;1;1;1;1;0];
    for ii = 1:groupsize
    rand_idx25(ii,1) = randsample(25,1,true, W); %,true,np/2 + b));
    W(rand_idx25(ii,1))=0;
    end
    
    %rand_idx25 = sort(rand_idx25)
    
    
    for ii = 1 : length(rand_idx25)
    rand_idx22(ii) = twenty2to25(rand_idx25(ii),'522');
    end
   
    
    e_idx{i} = zeros(size(selection_vector));
    e_idx{i}(rand_idx25) = 1;
end
%%
    V_CELL_sim = cell(6,1);
    V_CELL = cell(6,1);
    S_CELL = cell(6,1);
    
parfor idx = 1:length(V_CELL)-1
        [V_CELL{idx}, ~, S_CELL{idx}] = DFN_sim_casadi(p,...
            exp_num{1},Current_exp{1}(10:15), Time_exp{1}(10:15), ...
            Voltage_exp{1}(10:15), T_amb{1}, e_idx{idx}, theta, 1,Rc{1});

% see what happens if sensflag is 0 and we change senselect
    [V_CELL_sim{idx}] = DFN_sim_casadi(p,...
            exp_num{1},Current_exp{1}(10:15), Time_exp{1}(10:15), ...
            Voltage_exp{1}(10:15), T_amb{1}, e_idx{idx}, theta, 0,Rc{1});
end
%     
% idx = length(V_CELL);
%        [V_CELL{idx}] = DFN_sim_casadi(p,...
%             exp_num{1},Current_exp{1}(10:15), Time_exp{1}(10:15), ...
%             Voltage_exp{1}(10:15), T_amb{1}, 0, Nominal_param, 0,Rc{1});
%                [V_CELL_sim{idx}] = DFN_sim_casadi(p,...
%             exp_num{1},Current_exp{1}(10:15), Time_exp{1}(10:15), ...
%             Voltage_exp{1}(10:15), T_amb{1}, 0, Nominal_param, 0,Rc{1});
%%
figure(1)
hold on
for i = 1 : length(V_CELL)
    plot(V_CELL{i},'LineWidth',2)
    pause(1)
end
%%

disp('%%%%%%%%%%%%%%%')
for j = 1:length(V_CELL)
V_CELL_sim{j}' == V_CELL{j}'
disp('%%%%%%%%%%%%%%%')
end
%%
clear Vcomp
i = 1;
Vcomp = V_CELL{i};
disp('%%%%%%%%%%%%%%%')
for j = 1:length(V_CELL)
    j
Vcomp' == V_CELL{j}'
disp('%%%%%%%%%%%%%%%')
end
%V_CELL{1} == V_CELL{6}

%%
indvec = [1:4,7:19,21:25]
run param\params_nominal
viscomp = [cell2mat(e_idx),theta == Nominal_param];
visred = viscomp(indvec,:)
sum(and(visred(:,1:5) == 1, 1==visred(:,6)))


%%
indvec = [1:4,7:19,21:25]


L = length(indvec)
    V_CELL = cell(L+1,1);
    S_CELL = cell(L+1,1);
parfor idx = 1:L
    e_idx = zeros(size(selection_vector));
    e_idx(indvec(i)) = 1;
    
%         [V_CELL{idx}, ~, S_CELL{idx}] = DFN_sim_casadi(p,...
%             exp_num{1},Current_exp{1}(10:12), Time_exp{1}(10:12), ...
%             Voltage_exp{1}(10:12), T_amb{1}, e_idx, theta, 0,Rc{1});
            
        [V_CELL{idx}] = DFN_sim_casadi(p,...
            exp_num{1},Current_exp{1}(10:12), Time_exp{1}(10:12), ...
            Voltage_exp{1}(10:12), T_amb{1}, e_idx, theta, 0,Rc{1});
end

idx = L+1;
       [V_CELL{idx}] = DFN_sim_casadi(p,...
            exp_num{1},Current_exp{1}(10:12), Time_exp{1}(10:12), ...
            Voltage_exp{1}(10:12), T_amb{1}, 0, Nominal_param, 0,Rc{1});
%%

