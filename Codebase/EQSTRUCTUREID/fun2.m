function [ F ] = fun2(x,averageocp,N,Up,Un,p_soc,n_soc)
%THIS IS THE OBJECTIVE FUNCTION (IT RETURNS THE VALUE YOU WANT TO MINIMIZE
%IT IS SET UP TO RETURN THE MAXIMUM ERROR. IT CAN BE CHANGED TO RETURN THE 
%RMSE IF THAT IS WHAT YOU WANT TO OPTIMIZE... OR IT CAN BE CHANGED TO
%ANYTHING THAT YOU WANT TO OPTIMIZE. 

%NOTE: This function is called a lot so the shorter the runtime the better
%Hence I use a next-neighbor appoximation in the interp function. 

%Also it takes many values as inputs so that you don't have to re-load
%data. This same goal could be accomplished with global variables.

%% Set up SOC (linearly spaced-see pdf on eq struct. param. ID)
thetaN=linspace(x(1),x(2),N); %SOC of anode 
thetaP=linspace(x(3),x(4),N); %SOC of Cathode

%% Interpolates OCP
Upp = interp1(p_soc,Up,thetaP); %
Unn = interp1(n_soc,Un,thetaN); %

%% Find OCV from eq. model
ocv=Upp - Unn; 

%% Find error of eq model from real data for ocp
nnn=round(N/5);
range=[1:173129];
%range=[9115:N-nnn];
ocpe=averageocp(range)-ocv(range)'; %calculate error at each point
e2=ocpe.^2;
rmse=sqrt(sum(e2)/173129);
e=averageocp-ocv';
%% SOC error
%soc of
SOCexp=[0:1/(N-1):1];
SOCextrap=[1+1/(N-1):1/(N-1):1.005];
Extrappoly=polyfit(SOCexp(182260:end),averageocp(182260:end)',1);
OCPextrap=polyval(Extrappoly,SOCextrap);
[x, index] = unique([averageocp',OCPextrap]); 
y=[SOCexp,SOCextrap];
SOCmodel=interp1(x, y(index),ocv(173129:end));
SOCe=SOCmodel-SOCexp(173129:end);
FracErr=sum(ge(abs(ocpe),.005)/N);
F=5*1*FracErr+0*1*max(abs(ocpe))+0*max(abs(SOCe))+0*rmse;% Returns Max of error
end

