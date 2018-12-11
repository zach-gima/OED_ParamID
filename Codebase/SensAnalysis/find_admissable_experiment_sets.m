function [STSnorm,Sens_Mag,sens,NT_mat,exp_ind,num_exp] = find_admissable_experiment_sets(p,params,removed,Np,path,num_exp,r,bounds,output_folder)

    %% Normalization
    normalize_sens_factor = normalizesens(p,bounds);

    %%
    A=zeros(Np,num_exp);
    Sens_orth=zeros(Np,num_exp);
    Sens_Mag=zeros(Np,num_exp);
    AllR=cell(1,num_exp);
    AllQ=cell(1,num_exp);
    AllE=cell(1,num_exp);
    AllI=cell(1,num_exp);
    STS=cell(1,num_exp);
    STSnorm=cell(1,num_exp);
    sens = cell(1,num_exp);
    
    % experiment params
    T=zeros(1,num_exp);
    V_0_mat=zeros(1,num_exp);
    Sign=zeros(1,num_exp);
    NT_mat=zeros(1,num_exp);
    exp_ind = zeros(1,num_exp);
    

    mu = 10^-4;
    for i=1:num_exp
        t=r(i).name;

        if and(min(t(end-3:end) == '.mat'),t(1)~='e')
%             disp(t);
            index=str2double(t(1:end-4));
            exp_ind(i) = index; % use for storing experiment numbers to find later
%             disp(index);
            [T(i),V_0_mat(i),Sign(i),NT_mat(i),A(:,i) ,Sens_Mag(:,i),...
                Sens_orth(:,i),AllR{i},AllQ{i},AllE{i}, AllI{i},...
                STS{i},STSnorm{i}, sens{i}]=...
                extractexpinfo(t,path, normalize_sens_factor, Np, mu);
        else
            continue
        end

    end

    %% Plot Cardinality
    fig = plotcardinality(A, params,removed,Np);
    savefig(fig,strcat(output_folder,'cardinality_initial_Tmax60'));
    saveas(fig,strcat(output_folder,'cardinality_initial_Tmax60.png'));

    %% Plot Orth Sensitivity
    fig = plotcsensorth(Sens_orth, params, removed, -15, Np);
    savefig(fig,strcat(output_folder,'orthsens_initial_Tmax60'));
    saveas(fig,strcat(output_folder,'orthsens_initial_Tmax60.png'))

    %% Find Subsets
    %subsets=findsubsets(A, params)
    %%
%     save('/Users/ztakeo/Documents/GitHub/OED_ParamID/SensResults/exp_info_Tmax45.mat','A','AllE','AllQ','AllR','Sens_Mag','Sens_orth','T',...
%          'V_0_mat', 'Sign', 'NT_mat','AllI','sens','exp_ind','STSnorm','num_exp');
    save('/Users/ztakeo/Documents/GitHub/OED_ParamID/SensResults/exp_info_Tmax60.mat','A','AllE','AllQ','AllR','Sens_Mag','Sens_orth','T',...
         'V_0_mat', 'Sign', 'NT_mat','AllI','sens','exp_ind','STSnorm','num_exp');
end