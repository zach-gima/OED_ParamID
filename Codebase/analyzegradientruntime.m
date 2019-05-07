run param/params_nominal
init_cond = 'nom';
theta_0 = Nominal_param;

run param/params_NCA % loads p struct
run param/params_bounds % loads bounds struct (has fields bounds.min and bounds.max)
run param/params_truth % loads truth_param array

selection_vector = [1;1;1;1;0;0;1;1;1;1;1;1;1;1;1;1;1;1;1;0;1;1;1;1;1];
sel_k = find(selection_vector);

input_folder = strcat('InputLibrary/Experimental/');
filename_input = strcat(input_folder,'V_sim_G2G1.mat');
Inputs = load(filename_input);
Current_exp = Inputs.Current_exp;
Time_exp = Inputs.Time_exp; %% room for memory improvement
Voltage_exp = Inputs.V_LM_CELL;
T_amb = Inputs.T_amb_sim; % note, comes in celcius
exp_num = Inputs.exp_num;
Rc = Inputs.Rc;

time_S_save = [];
e_idx = zeros(1,25);

for i=sel_k'
e_idx(i) = 1
tic
for idx = 1
[V_CELL,~,Sens] = DFN_sim_casadi(p,...
            exp_num{idx},Current_exp{idx}(1:1), Time_exp{idx}(1:1), ...
            Voltage_exp{idx}(1:1), T_amb{idx}, e_idx, theta_0, 1,Rc{idx});
end
time_save = [time_save,toc];
hold on
plot(time_save);
drawnow
end

maxs = max(time_save);
mins = min(time_save);

time_normalized= (time_save - mins)/(maxs-mins);