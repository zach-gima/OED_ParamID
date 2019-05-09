    clear e_idx rand_idx25
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
    V_CELL = cell(5,1);
    S_CELL = cell(5,1);
    
parfor idx = 1:5
        [V_CELL{idx}, ~, S_CELL{idx}] = DFN_sim_casadi(p,...
            exp_num{1},Current_exp{1}(1:end), Time_exp{1}(1:end), ...
            Voltage_exp{1}(1:end), T_amb{1}, e_idx{idx}, theta, 1,Rc{1});
end
    
%%
figure
hold on
for i = 1 : 5
    plot(V_CELL{i})
end