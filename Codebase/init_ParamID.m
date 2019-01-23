%% init_paramID function: used to initialize variables for Param_ID function
% Post collinearity clustering & sensitivity analysis, only identifying
% 14 params:

% By Zach Gima 2018-11-9

%% Parameters [ZTG Change]

    %   UNCERTAIN PARAMETERS, theta
    %   G1: 2; G2: 6, G3: 5, G4: 5, Eq: 3 (Year 2)

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
    
function [filename_input_vector,filename_output_vector,selection_vector,ci_select,ci_input_vector] = init_ParamID(baseline,init_cond,input_folder,output_folder)
    
    %% Set Inputs, Parameters to Identify, and Inputs for Calculating Confidence Intervals.
    % These values are set differently based on 

    % Input Filenames (Year 2, post-collinearity and noise threshold clustering/elimination) 
    % Starting off just trying 2 groups of params
    filename_input_vector = cell(2,1);
    filename_input_vector{1} = strcat(input_folder,'V_sim_G1.mat');
    filename_input_vector{2} = strcat(input_folder,'V_sim_G2G1.mat');
    
    % Create output filenames and selection vector depending on approach
    filename_output_vector = cell(2,1);
    selection_vector = zeros(25,2);
    
    ci_select = cell(2,1);
    ci_input_vector = cell(2,1);

    % Baseline A: Full Parameter Set (1 Group)
    if strcmp(baseline{1},'full_') == 1 
        %Set output filename
        filename_output_vector{end} = strcat(output_folder,baseline{1},init_cond,'.mat');
        
        %Set selection vector
        selection_vector(:,end) = [1;1;1;1;0;0;1;1;1;1;1;1;1;1;1;1;1;1;1;0;1;1;1;1;1];
        
        ci_select{end} = find(selection_vector(:,end));
        ci_input_vector{end} = strcat(input_folder,'V_sim_G2G1.mat');
      
    % Baseline B: Collinearity Only (1 Group)
    elseif strcmp(baseline{1},'collinearity_') == 1
        filename_output_vector{end} = strcat(output_folder,baseline{1},init_cond,'.mat');

        %Set selection vector
        selection_vector(:,end) = [1;1;1;1;0;0;1;1;1;0;1;1;0;1;1;0;0;1;0;0;1;1;1;1;0];
        
        ci_select{end} = find(selection_vector(:,end));
        ci_input_vector{end} = strcat(input_folder,'V_sim_G2G1.mat');
        
    elseif strcmp(baseline{1},'OED_') == 1 % Baseline C: Collinearity + Sensitivity (2 Groups, can write for 1 as well)
        % Uncomment grouping approach you want to use
        %%%%%%%%%%%%%%%%%% All at once %%%%%%%%%%%%%%%%%% 
%         filename_output_vector{2} = strcat(output_folder,baseline{2},init_cond,'.mat');
% 
%         % Selection Vector
%         selection_vector(:,2) = [1;1;1;1;0;0;0;0;1;0;1;1;0;0;1;0;0;1;0;0;1;1;1;1;0]; % all together
%             
%         ci_select{2} = selection_vector(:,2);
%         ci_input_vector{2} = strcat(input_folder,'V_sim_G2G1.mat');
       
        %%%%%%%%%%%%%%%%%% Cumulative %%%%%%%%%%%%%%%%%% 
        %Set output filename
        filename_output_vector{1} = strcat(output_folder,baseline{1},init_cond,'.mat');
        filename_output_vector{2} = strcat(output_folder,baseline{1},init_cond,'.mat');
        
        % Selection vector (Year 2, post-collinearity and noise threshold clustering/elimination) 
        % Starting off just trying 2 groups of params
        % [ZTG Updated 2018-11-27]
        selection_vector(:,1) = [1;1;1;1;0;0;0;0;1;0;0;0;0;0;1;0;0;0;0;0;0;0;0;0;0]; %G1
        selection_vector(:,2) = [1;1;1;1;0;0;0;0;1;0;1;1;0;0;1;0;0;1;0;0;1;1;1;1;0]; %G2
        
        ci_select{1} = find(selection_vector(:,1));
        ci_select{2} = find(selection_vector(:,2) - selection_vector(:,1));

        ci_input_vector{1} = strcat(input_folder,'V_sim_G1.mat');
        ci_input_vector{2} = strcat(input_folder,'V_sim_G2.mat');
        
    else
        error('Incorrect baseline defined. Please check variable "baseline" ');
    end

    %% Confidence Interval variables
    % Create vector of indices for calculating confidence intervals
    % Will be calculating confidence intervals separately for each group
%     ci_select = cell(4,1);
%     ci_select{1} = find(selection_vector(:,1));
%     ci_select{2} = find(selection_vector(:,2) - selection_vector(:,1));
%     ci_select{3} = find(selection_vector(:,3) - selection_vector(:,2));
%     ci_select{4} = find(selection_vector(:,4) - selection_vector(:,3));
    
    % CI Input filenames
    % Note: for calculating confidence intervals, will only use the oed-cvx
    % selected inputs for that respective group
%     ci_input_vector = cell(4,1);
%     ci_input_vector{1} = strcat(input_folder,'V_sim_G1.mat');
%     ci_input_vector{2} = strcat(input_folder,'V_sim_G2.mat');
%     ci_input_vector{3} = strcat(input_folder,'V_sim_G3.mat');
%     ci_input_vector{4} = strcat(input_folder,'V_sim_G4.mat');
end