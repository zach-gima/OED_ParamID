fs = 25;

error = [129; 86.5]; %mV
time = [0; 45]; % hours

figure('Position', [100 100 900 700])
plot(time,error,'-o','LineWidth',3,'MarkerSize',15,'MarkerEdgeColor','r');
title('SQP Voltage RMSE');
xlabel('Wall Clock Time (Hr)')
ylabel('Voltage (mV)')
set(gca,'Fontsize',fs)
