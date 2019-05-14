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
%%
deltathetanorm = Param_norm_save(:,2:end) - Param_norm_save(:,1:end-1);
plot(sum(abs(deltathetanorm)))
figure
plot(Cost_save)
%%
ind = size(Param_norm_save);
figure
for i = 1:ind(2)-1
    hold off
plot(Param_norm_save(:,i),[1:22],'MarkerFaceColor','b','MarkerSize',15,'Marker','square','LineStyle','none')
hold on
plot(Param_norm_save(:,i)+deltathetanorm(:,i),[1:22],'MarkerFaceColor','r','MarkerSize',7,'Marker','square','LineStyle','none')
input('hi')
end


%% Plot Voltage
figure
plot(Voltage_truth_save,'k','LineWidth',2)
hold on
plot(Voltage_save(:,end))

% Cost_save(9)
% ans =
% 17.4315
% Cost_save(10)
% ans =
% 18.0990