%% init_paramID function: used to initialize variables for Param_ID function
% By Zach Gima 2018-3-19

%% Parameters [ZTG Change]

    %   UNCERTAIN PARAMETERS, theta
    %   G1: 2; G2: 6, G3: 5, G4: 5, Eq: 3 (Year 1)

    %   (G2) 1  : D_s_n       => p.D_s_n0
    %   (G2) 2  : D_s_p       => p.D_s_p0
    %   (G1) 3  : R_s_n       => p.R_s_n
    %   (G1) 4  : R_s_p       => p.R_s_p
    %   (Eq) 5  : epsilon_s_n => p.epsilon_s_n
    %   (Eq) 6  : epsilon_s_p => p.epsilon_s_p
    %   (G4) 7  : sig_n       => p.sig_n
    %   (G4) 8  : sig_p       => p.sig_p
    %   (G2) 9  : D_e         => p.ElecFactorD --> multiply factor
    %   (G2) 10 : epsilon_e_n => p.epsilon_e_n
    %   (G4) 11 : epsilon_e_s => p.epsilon_e_s
    %   (G3) 12 : epsilon_e_p => p.epsilon_e_p
    %   (G2) 13 : kappa       => p.ElecFactorK --> multiply factor
    %   (G4) 14 : t_plus      => p.t_plus
    %   (G2) 15 : dactivity   => p.ElecFactorDA --> multiply factor
    %   (G3) 16 : k_n0        => p.k_n0
    %   (G4) 17 : k_p0        => p.k_p0
    %   (G3) 18 : R_f_n       => p.R_f_n
    %   (G3) 19 : R_f_p       => p.R_f_p
    %   (Eq) 20 : n_Li_s      => p.n_Li_s
    %   (G3) 21 : c_e0        => p.c_e
    %   () 22 : E.Dsn        => p.E.Dsn
    %   () 23 : E.Dsp        => p.E.Dsp
    %   () 24 : E.kn        => p.E.kn
    %   () 25 : E.kp        => p.E.kp
    
function [filename_input_vector,filename_output_vector,selection_vector,ci_select,ci_input_vector] = init_ParamID(approach,init_cond,input_folder,output_folder)
    
    %% ParamID variables
%     % Selection vector (Year 1)
%     selection_vector = zeros(21,4); %21 parameters
%     selection_vector(:,1) = [0;0;1;1;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0]; %G1
%     selection_vector(:,2) = [1;1;1;1;0;0;0;0;1;1;0;0;1;0;1;0;0;0;0;0;0]; %G2
%     selection_vector(:,3) = [1;1;1;1;0;0;0;0;1;1;0;1;1;0;1;1;0;1;1;0;1]; %G3
%     selection_vector(:,4) = [1;1;1;1;0;0;1;1;1;1;1;1;1;1;1;1;1;1;1;0;1]; %G4

    % Selection vector (Year 1)
    selection_vector = zeros(25,4); %25 parameters
    selection_vector(:,1) = [0;0;1;1;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0]; %G1
    selection_vector(:,2) = [1;1;1;1;0;0;0;0;1;1;0;0;1;0;1;0;0;0;0;0;0;0;0;0;0]; %G2
    selection_vector(:,3) = [1;1;1;1;0;0;0;0;1;1;0;1;1;0;1;1;0;1;1;0;1;0;0;0;0]; %G3
    selection_vector(:,4) = [1;1;1;1;0;0;1;1;1;1;1;1;1;1;1;1;1;1;1;0;1;0;0;0;0]; %G4
    
    % File I/O
    % Input filenames
    filename_input_vector = cell(4,1);
    filename_input_vector{1} = strcat(input_folder,'V_sim_G1.mat');
    filename_input_vector{2} = strcat(input_folder,'V_sim_G2G1.mat');
    filename_input_vector{3} = strcat(input_folder,'V_sim_G3G2G1.mat');
    filename_input_vector{4} = strcat(input_folder,'V_sim_G4G3G2G1.mat');

    % Create output filenames
    filename_output_vector = cell(4,1);

    if strcmp(approach{4},'all_') == 1
        filename_output_vector{4} = strcat(output_folder,approach{4},init_cond,'.mat');
    else % cumulative approach
        filename_output_vector{1} = strcat(output_folder,approach{1},init_cond,'.mat');
        filename_output_vector{2} = strcat(output_folder,approach{2},init_cond,'.mat');
        filename_output_vector{3} = strcat(output_folder,approach{3},init_cond,'.mat');
        filename_output_vector{4} = strcat(output_folder,approach{4},init_cond,'.mat');
    end

    %% Confidence Interval variables
    % Create vector of indices for calculating confidence intervals
    % Will be calculating confidence intervals separately for each group
    ci_select = cell(4,1);
    ci_select{1} = find(selection_vector(:,1));
    ci_select{2} = find(selection_vector(:,2) - selection_vector(:,1));
    ci_select{3} = find(selection_vector(:,3) - selection_vector(:,2));
    ci_select{4} = find(selection_vector(:,4) - selection_vector(:,3));
    
    % CI Input filenames
    % Note: for calculating confidence intervals, will only use the oed-cvx
    % selected inputs for that respective group
    ci_input_vector = cell(4,1);
    ci_input_vector{1} = strcat(input_folder,'V_sim_G1.mat');
    ci_input_vector{2} = strcat(input_folder,'V_sim_G2.mat');
    ci_input_vector{3} = strcat(input_folder,'V_sim_G3.mat');
    ci_input_vector{4} = strcat(input_folder,'V_sim_G4.mat');
end