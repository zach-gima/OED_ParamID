%% Parameter Bounds
Np=[.0423,.1480]; %Lower and upper Bounds on N+
Nm=[.0763,.267];  %Lower and upper Bounds on N-
nli=[.001:.01:.6];

%% Load Data (CHANGE THIS)
 %load('criticalvars')
 load('UCB_ocps')
 load('C_50_OCP_data_final.mat')

%% Assign Variable shorter names
Un = UCB.anode.ocp; %OCP for anode
Up = UCB.cathode.ocp; %OCP for Cathode
n_soc = UCB.anode.soc; %corresopnding SOC for anode
p_soc = UCB.cathode.soc; %corresopnding SOC for cathode

%% TEST INFO
averageocp=OCPdata';     %<--- KEY LINE (set equal to ocv data, 
                               %properly truncated
N=length(averageocp)
dt = 1;                  % given sample time of your test 
                         % (ASSUMES EQUAL SAMPLING INTERVALS)
Fc = 96487;              % Faraday's constant [Coulombs/mol]
I = 0.058;               % CURRENT OF EXPERIMENT CURRENTLY for BOSCH C/50. [A]
nnn=round(N/5);
range=[1:164049];
%range=[9115:N-nnn];
%% BEGIN ASIGNING OPTIMIZATION VAR VALUES
th_m_z=[(N-1)*dt*I/Fc/Nm(2)+.001,1]; %FINDS UPPER AND LOWER BOUNDS ON Anode SOC
th_p_z=[.1461,1-(N-1)*dt*I/Fc/Np(2)-.01]; %FINDS UPPER AND LOWER BOUNDS ON CATHODE SOC
step=.01;  %CHANGE                          %STEP RESOLUTION FOR GUESSING


theta_minus_zero=[th_m_z(1):step:th_m_z(2)];
theta_plus_end=[.1461+(N-1)*dt*I/Fc/Np(2):step:1];

%solve for theta_plus
U_n_0=interp1(n_soc,Un,theta_minus_zero);
U_p_end=interp1(p_soc,Up,theta_plus_end);

U_n_end=-min(averageocp)+U_p_end;
U_p_0=max(averageocp)+U_n_0;

theta_plus_zero=interp1(Up,p_soc,U_p_0);
theta_minus_end=interp1(Un,n_soc,U_n_end);

%%
N1=length(theta_minus_zero);
N2=length(theta_minus_end);

amtOSteps=N1*N2; %amount of steps

%% INITIALIZE MATRICIES FOR DATA STORAGE

xvals=cell(amtOSteps,1);
fvals=ones(amtOSteps,4);
Flag=ones(amtOSteps,1);
i=1;
%% CONSTRAINT MATRICIES FOR OPTIMIZATION
A=[1 0 0 0;-1 0 0 0;0 0 1 0;0 0 -1 0;0 0 0 1;-1 1 0 0;0 -1 0 0;0 0 1 -1];
B=[1;-.4114;.2596;-.1461;1;-.4104;0.01;-.7404];
%OPTIONS = optimoptions('fmincon','Algorithm','interior-point')
%OPTIONS = optimoptions('fmincon','Algorithm','sqp')
%OPTIONS = optimoptions('fmincon','Algorithm','trust-region-reflective')
%OPTIONS = optimoptions('fmincon','Algorithm','active-set')
%% ITERATE OVER INITIAL GUESSES (IMPROVEMENT On THIS SECTION IN THE WORKS)
j=1
count=0
ind=1
for i=1:N1
    for j=1:N2
        %interior point, SQP, active set, and trust region reflective
        xguess=[theta_minus_zero(i),theta_minus_end(j),theta_plus_zero(i),theta_plus_end(j)];
        [fvals(ind,1),fvals(ind,2),fvals(ind,3),fvals(ind,4)]=funiter(xguess,averageocp,N,Up,Un,p_soc,n_soc,range);
        %[ FracErr, MaxErr,SOCerr,rmse ]
        %[xvals{ind},fvals(ind),Flag(ind)]=fmincon(@(x) fun2(x,averageocp,N,Up,Un,p_soc,n_soc),xguess,A,B,[],[],[],[],[],OPTIONS);
        xvals{ind}=xguess;
         
        ind=ind+1
    end
    j=j+1; %REMOVE SEMI COLON TO GET A SENSe OF THE RUNTIME
end



%% FIND BEST OPTIMIZATION PARAMETERS AND MAX MIN ERROR
%f=find(fvals+1000==min(fvals+eq(Flag,2)*1000+eq(Flag,1)*1000+eq(Flag,-2)*100000+eq(Flag,0)*100000));
%[ FracErr, MaxErr,SOCerr,rmse ]
fvalflag=1
FracErr=fvals(:,1);
MaxErr=fvals(:,2);
SOCerr=fvals(:,3);
rmse=fvals(:,4);

%%
f=find(fvals(:,fvalflag)==min(fvals(:,fvalflag)))
fopt=fvals(f,fvalflag);
xopt=xvals{f};
x=xopt;


%%
thetaN=linspace(x(1),x(2),N); %SOC of anode 
thetaP=linspace(x(3),x(4),N); %SOC of Cathode

%% Interpolates OCP
Upp = interp1(p_soc,Up,thetaP); %
Unn = interp1(n_soc,Un,thetaN); %

%% Find OCV from eq. model
ocv=Upp - Unn; 

%% Find error of eq model from real data for ocp

ocper=averageocp(range)-ocv(range)'; %calculate error at each point
e2=ocper.^2;
rmse=sum(e2)/173129;
er=averageocp-ocv';
%% SOC error
%soc of
SOCexp=[0:1/(N-1):1];
SOCextrap=[1+1/(N-1):1/(N-1):1.005];
Extrappoly=polyfit(SOCexp(182260:end),averageocp(182260:end)',1);
OCPextrap=polyval(Extrappoly,SOCextrap);
[xocp, index] = unique([averageocp',OCPextrap]); 
y=[SOCexp,SOCextrap];
SOCmodel=interp1(xocp, y(index),ocv(range(end):end));
SOCer=SOCmodel-SOCexp(range(end):end);
e=[ocper',SOCer];

SOCe=max(abs(SOCer))
ocpe=max(abs(ocper))
%%
%Calculate N+,N-,n_li
disp('eq parameters')
Np=N*dt*I/Fc/(x(4)-x(3));
Nm=N*dt*I/Fc/(x(1)-x(2));
nli=Nm*x(1)+Np*x(3);
disp({'Np=',Np})
disp({'Nm=',Nm})
disp({'nli=',nli})

%%
disp('optimal anode and cathode SOC Ranges')
disp({'theta.anode(0)',x(1)})
disp({'theta.anode(end)',x(2)})
disp({'theta.cathode(0)',x(3)})
disp({'theta.cathode(end)',x(4)})
disp('objective value')
disp('minimum of max error')
disp(fopt)
x=xopt;
thetaN=linspace(x(1),x(2),N);
thetaP=linspace(x(3),x(4),N);
Upp = interp1(p_soc,Up,thetaP);
Unn = interp1(n_soc,Un,thetaN);
ocv=Upp - Unn;
figure
SOC=linspace(0,1,N);
plot(SOC,ocv)
hold on
plot(SOC,averageocp)
plot(SOC,averageocp+.005,'k-')
plot(SOC,averageocp-.005,'k-')

  figure
  hold on
    plot(n_soc,Un)
 plot([x(1),x(1)],[0,1.4],'k-')
 plot([x(2),x(2)],[0,1.4],'k-')
 set(gca,'FontSize', 16);
 title('anode')
 xlabel('SOC')
 ylabel('OCP')
 
   figure
  hold on
    plot(p_soc,Up)
 plot([x(3),x(3)],[3,4.5],'k-')
 plot([x(4),x(4)],[3,4.5],'k-')
 set(gca,'FontSize', 16);
  title('cathode')
 xlabel('SOC')
 ylabel('OCP')
 
%% cdf

off=ge(abs(ocper),.005);
reg=[0:.001:.1];
cdf=zeros(1,length(reg));
count=1;
for i=reg;
    cdf(count)=sum(le(abs(ocper),i));
    count=count+1;
end
figure
plot(cdf/length(er)*100,reg*1000);
hold on
plot([0 100],[5 5]);
plot([50 50],[0 max(reg*1000)]);
ylabel('Absoulute error mV')
xlabel('percentile')
title('CDF For Minimizing error > 5mV')
legend('CDF','5mV','50 percentile')



figure
plot(SOC(range),abs(ocper));
ylabel('Error')
xlabel('SOC')
title('Abs(Error) vs SOC for Minimizing error > 5mV')

%% cdf

off=ge(abs(ocper),.005);
reg=[0:.0001:.02];
cdf=zeros(1,length(reg));
count=1;
for i=reg;
    cdf(count)=sum(le(abs(SOCer),i));
    count=count+1;
end
figure
plot(cdf/length(er)*1000,reg*100);
hold on
%plot([0 100],[5 5]);
%plot([50 50],[0 max(reg*1000)]);