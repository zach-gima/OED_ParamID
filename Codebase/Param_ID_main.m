%% OED_ParamID fmincon implementation
% By: Zach Gima 2019-5-6
clc
clearvars
close all

datetime_initial = datetime('now','TimeZone','America/Los_Angeles');

%% instantiate global variables that will track iterations within fmincon
global history;
global searchdir;
history.x = [];
history.fval = [];
searchdir = [];

%% User Input 

% Load Nominal Parameter Set, Bounds, & True params
run param/params_NCA % loads p struct
run param/params_bounds % loads bounds struct (has fields bounds.min and bounds.max)

%%%%%%%%%%%%%%%   ParamID baseline (Uncomment to select)   %%%%%%%%%%%%%%%

% Baseline C: Collinearity + Sensitivity (2 Groups)
baseline = {'OED_'};
num_groups = 1; % Number of parameter groups

%%%%%%%%%%%%%%%   Parameter Initial Conditions (Uncomment to select)   %%%%%%%%%%%%%%%
run param/params_nominal
init_cond = 'nom';
theta_0 = Nominal_param;

%%%%%%%%%%%%%%%  File I/O (Set once)   %%%%%%%%%%%%%%%
% Input subfolder
input_folder = strcat('InputLibrary/Experimental/');

% Output subfolder
date_txt = strrep(datestr(datetime_initial), ':', '_');
% output_folder = strcat('/Users/ztakeo/Documents/GitHub/OED_ParamID/ID_results/EEC227C/',date_txt,'/');
% output_folder = strcat('C:/Users/Zach/Box Sync/HPC/HPC1/',date_txt,'/'); %HPC-1 Path
% output_folder = strcat('C:/Users/zgima/Box Sync/HPC/HPC2/',date_txt,'/'); %HPC-2 Path
% output_folder = strcat('/global/home/users/ztakeo/output/',date_txt,'/'); %Savio Path

mkdir(output_folder); %create new subfolder with current date in output_folder

%%% init_ParamID: initialize background stuff (variables, file i/o etc) based on the ParamID baseline and I.C.'s 
% [filename_input_vector,filename_output_vector,selection_vector,ci_select,ci_input_vector] = init_ParamID(baseline,init_cond,num_groups,input_folder,output_folder);
filename_input_vector{1} = strcat(input_folder,'V_sim_G2G1.mat');
filename_output_vector{1} = strcat(output_folder,baseline{1},'G2G1_',init_cond,'.mat');
% selection_vector = [1;1;1;1;0;0;0;0;1;0;1;1;0;0;1;0;0;1;0;0;1;1;1;1;0]; % 13 params, G2
% selection_vector = [1;1;1;1;0;0;0;0;1;0;0;0;0;0;1;0;0;0;0;0;0;0;0;0;0]; % 6 params, G1
% selection_vector = [1;1;1;1;0;0;1;1;1;1;1;1;1;1;1;1;1;1;1;0;1;1;1;1;1]; % 22 params
% selection_vector = [1;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0]; % debug

%% Display Simulation Info

%try/catch structure used to send email alert if program exits w/ error
% For saving errors:
error_filename = strcat(output_folder,'sim_log.txt');
diary(error_filename)

datetime_initial
fprintf('Initial Conditions: %s \n',init_cond);
fprintf('Baseline: %s \n',baseline{1});
fprintf('Number of Groups: %i \n',num_groups);
fprintf('Number of Parameters: %i \n',sum(selection_vector));

%% Call ParamID function
try
    for jj = 1:num_groups
        % Load Group Specific Inputs
        filename_input = filename_input_vector{jj};
        filename_output = filename_output_vector{jj};
        Inputs = load(filename_input); %Current, Voltage, Time, T_amb     
        sel_k = find(selection_vector);

        % Setup vectors related to parameter sensitivity selection
        selected_params = theta_0(sel_k);  % Goes from 21x1 vector to 18x1 (in all params selected scenario)
        num_param = size(selected_params,1); % number of param  

        % run fmincon     
        tic

        lb = bounds.min(sel_k);
        ub = bounds.max(sel_k);

        opt = optimoptions('fmincon','Display','iter','Algorithm','sqp',...
        'SpecifyObjectiveGradient',true,'Diagnostics','on','OutputFcn',@outfun);  %'CheckGradients',true,'FiniteDifferenceType','central',

        % create anonymous function to pass in other parameters needed for V_obj
        fh = @(x)V_obj(x,selection_vector, Inputs, p);

        [theta_ID,FVAL,EXITFLAG,fmincon_output,LAMBDA,GRAD,HESSIAN]= ...
           fmincon(fh,selected_params,[],[],[],[],lb,ub,[],opt);

        t_wallclock = toc;
        datetime_final = datetime('now','TimeZone','local')
        fprintf('Parameter ID complete. Finished in %i seconds. \n',t_wallclock);

        % Save data & send email
        save(filename_output,'theta_ID','fmincon_output','t_wallclock','history','searchdir');

        % matlabmail(recipient,subject,message,attachments)
        matlabmail('ztakeo@berkeley.edu','Parameter ID complete','',[]);
    end
catch e %e is an MException struct
    % An error will put you here.
    errorMessage = sprintf('%s',getReport( e, 'extended', 'hyperlinks', 'on' ))
    
%     % Save ParamID results even when erorr occurs
%     save(filename_output,'park0','paramID_out','LM_Iter','ci95_full');
    save(strcat('debug_',filename_output),'searchdir','history');
    diary off %stop logging command window
    
    matlabmail('ztakeo@berkeley.edu','Error encountered',errorMessage,error_filename);
end    

%% The purpose of this function is to save the iteration values 
%found and slightly modified from online documentation
function stop = outfun(x,optimValues,state)
    stop = false;
    global history;
    global searchdir;
    switch state
        case 'init'
            hold on
        case 'iter'
            % Concatenate current point and objective function
            % value with history. x must be a row vector.
            history.fval = [history.fval; optimValues.fval];
            history.x = [history.x; x'];
            % Concatenate current search direction with 
            % searchdir.
            searchdir = [searchdir;... 
                        optimValues.searchdirection'];
            plot(x(1),x(2),'o');
            % Label points with iteration number and add title.
            % Add .15 to x(1) to separate label from plotted 'o'
            text(x(1)+.15,x(2),... 
                num2str(optimValues.iteration));
            title('Sequence of Points Computed by fmincon');
        case 'done'
            hold off
        otherwise
    end
end
