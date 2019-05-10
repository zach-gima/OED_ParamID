%Debug Script
figure
plot(Cost_save)
%%
errorvec= Cost_save(2:end)-Cost_save(1:end-1)>0;

figure
plot(errorvec)

deltathetanorm = -Param_norm_save(:,1:end-1)+Param_norm_save(:,2:end)
changedparams = (deltathetanorm ~= 0);
didchang = sum(changedparams);

%%
deltatheta = -Param_save(:,1:end-1)+Param_save(:,2:end)
changedparams = (deltatheta ~= 0);
didchang = sum(changedparams);

%%
diff = (Voltage_save-Voltage_truth_save);
[sqrt(sum(diff.^2))',Cost_save]

%%
run param/params_NCA
run param/params_bounds
i = 1;
a=Param_save(:,i)<bounds.min;
b=Param_save(:,i)>bounds.max;
bounds.min<bounds.max

% Cost_save(9)
% ans =
% 17.4315
% Cost_save(10)
% ans =
% 18.0990