tic
clear
%% instantiate global variables that will track iterations within fmincon
global history;
global searchdir;
history.x = [];
history.fval = [];
searchdir = [];

%% get NCA parameters
run param/params_NCA
G1ind = [3,4]; %variable positions
G2ind = [1,2,7,8,11,13];  %variable positions
inds = G1ind; %variable positions for this case

%% make optimization options object
%this sets the algorithm, line search bounds
%Also sets an "output function" this was an addition for debugging that
%helps us look inside fmincon
opt = optimoptions('fmincon','Algorithm','active-set','RelLineSrchBnd',.5,'RelLineSrchBndDuration',7, ...
'Diagnostics','on','OutputFcn',@outfun); %sets 3 optimization options

%% get parameters for simulation
SensSelec = [1;1;1;1;0;0;1;1;1;1;1;1;1;1;1;1;1;1;1;0;1]; %for defining selk
sel_k = find(SensSelec); %vector for selecting 18 of 21 parameters
run param/params_truth
x0 = truth_param(sel_k); %initial guess
perturb = .9;
x0 = x0(inds) * perturb;

%% parse bounds
run param/params_bounds
lb = bounds.min(sel_k);%upper bound
ub = bounds.max(sel_k);%lower bound
lb = lb(inds);
ub = ub(inds);

%% Matlab documentations suggested to try this if the algorithm is taking weird steps
opt.TypicalX = sqrt(lb.*ub);

%% solve the problem 
% the inputs to fmincon are:
%the objective function, initial guess, bounds, and optimization object
%the only confusion might be the @g1fun handle. It maps the parameters theta 
%that are trying to be identified to the corresponding rmse. see @g1 fun
%for more
[X,FVAL,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN]= ...
    fmincon(@g1fun,x0,[],[],[],[],lb*.9999,ub*1.0001,[],opt)

totaltime = toc
%try running later with a shorter input




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