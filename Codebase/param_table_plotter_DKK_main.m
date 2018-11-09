%%% To run param_table_plotter independently from automated script, 
% (1) load final group Param ID results
% (2) copy/paste below and run:

function param_table_plotter_DKK_main(ci95_full,selection_vector,group,paramID_out,wait,theta_0_true,truth_param,bounds)
    % Indexing related variables
    n=length(paramID_out.save_param_org(1,:));
    names={'G1','G1,G2','G1,G2,G3','G1,G2,G3,G4'};
    for i = [1:n]
    %indicates which variables were identified (currently Group 4, so all but the eq. struct. params)
    sel_k = find(selection_vector(:,1)); %indices of the selected parameter values

    park0 = paramID_out.save_param_org(:,i); %loads the identified parameters. 

    theta_0 = zeros(21,1);
    theta_0(sel_k) = park0;
    drawnow nocallbacks
    param_table_animator_DKK(ci95_full,theta_0,theta_0_true,truth_param,bounds);
    title(names{group});
    drawnow;
    pause(wait);
    %@Zach Uncomment the line below to advance one step on eyboard input
    %input('advance');
    end
    %Zach note for Dylan: 
    % When I first thought of animating the function, I figured I would need to
    % call the param_table_plotter_ZTG in a for loop. Ignore if you find a
    % better way, but if that's the route you take then you'll also need to
    % update (1) selection_vector and (2) park0 within the for loop depending
    % on which Group of parameters you're animating. Or if you're animating all
    % of the parameters at the same time, then maybe it's not an issue.
end