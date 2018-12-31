%Animator
function animate_DKK(output_folder,selection_vector,truth_param)
    wait = 0;
        
    %In future load selection_vector as well after saving to .mat file
    flag = 1; % for gif making purposes
    load(strcat(output_folder,'G1_nom.mat'),'ci95_full','paramID_out','theta_0_true','bounds');
    group = 1;
    param_table_plotter_DKK_main(ci95_full,selection_vector,group,paramID_out,wait,theta_0_true,truth_param,bounds,flag,output_folder)

    flag = 0; % for gif making purposes
    load(strcat(output_folder,'G2G1_nom.mat'),'ci95_full','paramID_out','theta_0_true','bounds');
    group = 2;
    param_table_plotter_DKK_main(ci95_full,selection_vector,group,paramID_out,wait,theta_0_true,truth_param,bounds,flag,output_folder)

%     load(horzcat(output_folder,'G3G2G1_nom.mat'))
%     selection_vector = [1;1;1;1;0;0;0;0;1;1;0;1;1;0;1;1;0;1;1;0;1];
%     group = 3;
%     param_table_plotter_DKK_main(ci95_full,selection_vector,group,paramID_out,wait)
% 
%     load(horzcat(output_folder,'G4G3G2G1_nom.mat'))
%     selection_vector = [1;1;1;1;0;0;1;1;1;1;1;1;1;1;1;1;1;1;1;0;1];
%     group = 4;
%     param_table_plotter_DKK_main(ci95_full,selection_vector,group,paramID_out,wait)
end