%%% To run param_table_plotter independently from automated script, 
% (1) load final group Param ID results
% (2) copy/paste below and run:

function param_table_plotter_DKK_main(ci95_full,selection_vector,group,paramID_out,wait,theta_0_true,truth_param,bounds,flag,output_folder)
    % Indexing related variables
    num_iter=length(paramID_out.save_param_org(1,:));
    names={'G1','G1,G2','G1,G2,G3','G1,G2,G3,G4'};
    
    %indicates which variables were identified
    sel_k = find(selection_vector(:,group)); %indices of the selected parameter values
    
    % initialize Final param vector; used to update params
    Final_param = zeros(length(theta_0_true),1);
    
    %Setup .gif saving info
%     h = figure(2);
%     set(gcf,'Position',[100 100 900 700],'PaperPositionMode','auto');
% 
%     axis tight manual % this ensures that getframe() returns a consistent size
    filename = strcat(output_folder,'dartboard.gif');

    % Animate each L-M iteration of Param updates
    for i = 1:num_iter

        park0 = paramID_out.save_param_org(:,i); %loads the identified parameters. 
        Final_param(sel_k) = park0; %update Final_param with current iteration

        %     theta_0 = zeros(21,1);
        %     theta_0(sel_k) = park0;
        drawnow nocallbacks

        % Plot. Note param_table_plotter should always use final group's
        % selection_vector
        g = param_table_plotter_ZTG('animate',selection_vector(:,end),ci95_full,Final_param,theta_0_true,truth_param,bounds);

        title(names{group});
        drawnow;
        pause(wait);

        % Capture the plot as an image 
        frame = getframe(g); 
        im = frame2im(frame); 
        [imind,cm] = rgb2ind(im,256); 
        
        % Write to the GIF File 
      if i == 1 && flag == 1
          imwrite(imind,cm,filename,'gif','Loopcount',inf); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append'); 
      end 
        %@Zach Uncomment the line below to advance one step on eyboard input
        %         input('advance');
    end
end