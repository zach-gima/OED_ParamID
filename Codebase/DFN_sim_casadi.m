%%%%%%%%%% ZTG Edit of sens_lm_cal_noRc.m

%%%NOTE: DFN_casadi must take in a SelecParam vector of 18x1 (i.e.
%%%excluding the equilibrium parameters). Otherwise, park0 will be assigned
%%%incorrectly (parameters will be assigned the valeus for other
%%%parameters, and IDAS will throw an error because  the integration will
%%%yield impossible values

% Modified July 9, 2018 by Zach Gima (edit of SHP sens_lm_cal_noRc function).
% Descr: This function simulates the Doyle-Fuller-Newamn model using a CasADi
% automatic differentiation method. This modification merged functionality
% of DFN_sim_noRc + sens_lm_cal_noRc. Sensitivity of the output voltage
% w.r.t. parameters in the DFN can be calculated if the feature is turned
% on (SensFlag == 1). Otherwise (SensFlag == 0), only voltage is simulated

% T_amb expected in Celsius

function [v_sim,alg_states,varargout] = DFN_sim_casadi(p, Current_exp, Time_exp, Voltage_exp, T_amb, SensSelec, SelecParam, SensFlag,Rc) % [ZTG change] removed Rc for no model-to-model comparison

    addpath('/Users/ztakeo/Documents/MATLAB/casadi') % Mac Laptop
%     addpath('C:/Users/Zach/Documents/MATLAB/casadi_windows') % HPC-1
%     addpath('C:/Users/zgima/Documents/MATLAB/casadi_windows') % HPC-2
%     addpath('/global/home/users/ztakeo/modules/casadi-matlab');    % For Savio
    import casadi.*

    %% Load Electrochemical Model Parameters

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Battery Model : SPMe
    % disp('Battery Model Setting: DFN')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    str={'p.D_s_n0', 'p.D_s_p0', 'p.R_s_n', 'p.R_s_p', 'p.epsilon_s_n', 'p.epsilon_s_p', 'p.sig_n', 'p.sig_p' ...
         'p.ElecFactorD', 'p.epsilon_e_n', 'p.epsilon_e_s', 'p.epsilon_e_p', 'p.ElecFactorK', 'p.t_plus', ...
         'p.ElecFactorDA','p.k_n0','p.k_p0', 'p.R_f_n', 'p.R_f_p','p.n_Li_s', 'p.c_e0', 'p.E.Dsn','p.E.Dsp','p.E.kn','p.E.kp'};

    selection_vector = SensSelec;
    sel_k = find(selection_vector);
    np = length(sel_k);
    
    % assign
    D_s_n0 = p.D_s_n0;
    D_s_p0 = p.D_s_p0;
    R_s_n = p.R_s_n;
    R_s_p = p.R_s_p;
    epsilon_s_n = p.epsilon_s_n;
    epsilon_s_p = p.epsilon_s_p;
    sig_n = p.sig_n;
    sig_p = p.sig_p;
    ElecFactorD = p.ElecFactorD;
    epsilon_e_n = p.epsilon_e_n;
    epsilon_e_s = p.epsilon_e_s;
    epsilon_e_p = p.epsilon_e_p;
    ElecFactorK = p.ElecFactorK;
    t_plus = p.t_plus;
    ElecFactorDA = p.ElecFactorDA;
    k_n0 = p.k_n0;
    k_p0 = p.k_p0;
    R_f_n = p.R_f_n;
    R_f_p = p.R_f_p;
    n_Li_s = p.n_Li_s;
    c_e0 = p.c_e0;
    E.Dsn = p.E.Dsn;
    E.Dsp = p.E.Dsp;
    E.kn = p.E.kn;
    E.kp = p.E.kp;
    
    a_s_n = 3*p.epsilon_s_n / p.R_s_n;
    a_s_p = 3*p.epsilon_s_p / p.R_s_p;
    epsilon_f_n = 1 - p.epsilon_s_n - p.epsilon_e_n;
    epsilon_f_p = 1 - p.epsilon_s_p - p.epsilon_e_p;

    % Save these parameters into cell_var for comparison later
    cell_var = {D_s_n0,...
                D_s_p0,...
                R_s_n,...
                R_s_p,...
                epsilon_s_n,...
                epsilon_s_p,...
                sig_n,...
                sig_p,...
                ElecFactorD,...
                epsilon_e_n,...
                epsilon_e_s,...
                epsilon_e_p,...
                ElecFactorK,...
                t_plus,...
                ElecFactorDA,...
                k_n0,...
                k_p0,...
                R_f_n,...
                R_f_p,...
                n_Li_s,...
                c_e0,...
                E.Dsn,...
                E.Dsp,...
                E.kn,...
                E.kp,...
                a_s_n,...
                a_s_p,...
                epsilon_f_n,...
                epsilon_f_p};

    % Convert CasADi variables
    str_tmp = str;
    for i=1:length(sel_k)
        str_tmp{sel_k(i)}=SX.sym(str{sel_k(i)},1);
        cell_var{sel_k(i)} = str_tmp{sel_k(i)};
    end        

    % assign CasADi variables
    p.D_s_n0 = cell_var{1};
    p.D_s_p0 = cell_var{2};
    p.R_s_n = cell_var{3};
    p.R_s_p = cell_var{4};
    p.epsilon_s_n = cell_var{5};
    p.epsilon_s_p = cell_var{6};
    p.sig_n = cell_var{7};
    p.sig_p = cell_var{8};
    p.ElecFactorD = cell_var{9};
    p.epsilon_e_n = cell_var{10};
    p.epsilon_e_s = cell_var{11};
    p.epsilon_e_p = cell_var{12};
    p.ElecFactorK = cell_var{13};
    p.t_plus = cell_var{14};
    p.ElecFactorDA = cell_var{15};
    p.k_n0 = cell_var{16};
    p.k_p0 = cell_var{17};
    p.R_f_n = cell_var{18};
    p.R_f_p = cell_var{19};
    p.n_Li_s = cell_var{20};
    p.c_e0 = cell_var{21};
    p.E.Dsn = cell_var{22};
    p.E.Dsp = cell_var{23};
    p.E.kn = cell_var{24};
    p.E.kp = cell_var{25};

    p.a_s_n = 3*p.epsilon_s_n / p.R_s_n;
    p.a_s_p = 3*p.epsilon_s_p / p.R_s_p;
    p.epsilon_f_n = 1 - p.epsilon_s_n - p.epsilon_e_n;
    p.epsilon_f_p = 1 - p.epsilon_s_p - p.epsilon_e_p;

    %% Input charge/discharge Current Data %%
    % % Current | Positive <=> Discharge, Negative <=> Charge

    cn_low = 44.52; % Minimum stochiometry of Anode (CELL_SOC=0)
    cn_high = 34310.08; % Maximum stoichiometry of Anode (CELL_SOC=1)

    cp_low = 10118.4; % Minimum stochiometry of Cathode (CELL_SOC=1)
    cp_high = 46053; % Maximum stochiometry of Cathode (CELL_SOC=0)

    Delta_cn = cn_high - cn_low;
    Delta_cp = cp_high - cp_low;

    OneC = min(p.epsilon_s_n*p.L_n*Delta_cn*p.Faraday/3600, p.epsilon_s_p*p.L_p*Delta_cp*p.Faraday/3600);
    p.OneC = OneC;% [A/m^2]

    p.delta_t = 1;

    %%%%%%%%%%%%%%% Conference Test %%%%%%%%%%%%%%%%%%%%%%%%%
    % load('InputLibrary/1Ccrg_dchrg.mat')
    V0 = Voltage_exp(1);
    t = 1:length(Time_exp);
    I = -Current_exp/p.Area;
    p.R_c = Rc;
    NT = length(t);

    %% (DFN Code Copy) Initial Conditions & Preallocation

    % Given V0, obtain 
    [csn0,csp0] = init_cs(p,V0);
    V0 = refPotentialCathode(p,csp0/p.c_s_p_max)-refPotentialAnode(p,csn0/p.c_s_n_max);

    % Electrolyte concentration
    ce0 = p.c_e0;

    % Temperature
%     T0 = p.T_amb; % [ZTG Change] Taken in as an input from the input
%     profile in Year 2

    % Vector lengths
    Ncsn = p.PadeOrder * (p.Nxn-1);
    Ncsp = p.PadeOrder * (p.Nxp-1);
    Nn = p.Nxn - 1;
    Ns = p.Nxs - 1;
    Np = p.Nxp - 1;
    Nx = p.Nx - 3;

    c_s_n0 = zeros(p.PadeOrder,1);
    c_s_p0 = zeros(p.PadeOrder,1);

    %%%%% Initial condition based on Jordan form
    c_s_n0(3) = csn0;
    c_s_p0(3) = csp0;

    c_s_n = zeros(Ncsn,NT);
    c_s_p = zeros(Ncsp,NT);

    c_s_n(:,1) = repmat(c_s_n0, [Nn 1]);
    c_s_p(:,1) = repmat(c_s_p0, [Nn 1]);

    % Electrolyte concentration
    %c_e = zeros(Nx,NT); %CASADI CHANGE
    c_e = SX.zeros(Nx,NT);
    c_e(:,1) = ce0 * ones(Nx,1);

    % Temperature
    % T = zeros(NT,1);
    % T(1) = T0;
%     T10 = p.T_amb;
%     T20 = p.T_amb;
    T10 = T_amb(1) + 273.15; % T_amb expected in Celsius
    T20 = T_amb(1) + 273.15;

    p.T_amb = T_amb(1) + 273.15;
    
    % Solid Potential
    Uref_n0 = refPotentialAnode(p, csn0(1)*ones(Nn,1) / p.c_s_n_max);
    Uref_p0 = refPotentialCathode(p, csp0(1)*ones(Np,1) / p.c_s_p_max);

    phi_s_n = zeros(Nn,NT);
    phi_s_p = zeros(Np,NT);
    phi_s_n(:,1) = Uref_n0;
    phi_s_p(:,1) = Uref_p0;

    % Electrolyte Current
    i_en = zeros(Nn,NT);
    i_ep = zeros(Np,NT);

    % Electrolyte Potential
    phi_e = zeros(Nx+2,NT);

    % Molar Ionic Flux
    jn = zeros(Nn,NT);
    jp = zeros(Np,NT);
    
    % Volume average concentration
    c_avg_n = zeros(Nn,NT);
    c_avg_n(:,1) = repmat(csn0, [Nn 1]);

    % SOC (Bulk Anode SOC)
    SOC = zeros(NT,1);
    SOC(1) = (mean(c_avg_n(:,1)) - cn_low) / (cn_high - cn_low);%soc00;%(mean(c_avg_n(:,1)) - cn_low) / (cn_high - cn_low);

    % Initial Conditions
    x0_nom = [c_s_n(:,1); c_s_p(:,1); c_e(:,1); T10; T20];

    z0 = [phi_s_n(:,1); phi_s_p(:,1); i_en(:,1); i_ep(:,1);...
          phi_e(:,1); jn(:,1); jp(:,1)];

    %% (DFN Code Copy) Pre-computation
    % Electrolyte concentration matrices
    [M1n,M2n,M3n,M4n,M5n, M1s,M2s,M3s,M4s, M1p,M2p,M3p,M4p,M5p, C_ce] = c_e_mats(p);

    p.ce.M1n = M1n;
    p.ce.M2n = M2n;
    p.ce.M3n = M3n;
    p.ce.M4n = M4n;
    p.ce.M5n = M5n;

    p.ce.M1s = M1s;
    p.ce.M2s = M2s;
    p.ce.M3s = M3s;
    p.ce.M4s = M4s;

    p.ce.M1p = M1p;
    p.ce.M2p = M2p;
    p.ce.M3p = M3p;
    p.ce.M4p = M4p;
    p.ce.M5p = M5p;

    p.ce.C = C_ce;

    rM3 = [Nn; Ns; Np];
    cM3 = rM3';
    p.ce.M3 = sparse(blkdiagFast(rM3, cM3, p.ce.M3n, p.ce.M3s, p.ce.M3p));

    clear M1n M2n M3n M4n M5n M1s M2s M3s M4s M1p M2p M3p M4p M5p C_ce rM1 cM1;

    % Solid Potential
    [F1_psn,F1_psp,F2_psn,F2_psp,G_psn,G_psp,...
        C_psn,C_psp,D_psn,D_psp] = phi_s_mats(p);
    p.F1_psn = F1_psn;
    p.F1_psp = F1_psp;
    p.F2_psn = F2_psn;
    p.F2_psp = F2_psp;
    p.G_psn = G_psn;
    p.G_psp = G_psp;
    p.C_psn = C_psn;
    p.C_psp = C_psp;
    p.D_psn = D_psn;
    p.D_psp = D_psp;

    clear F1_psn F1_psp F2_psn F2_psp G_psn G_psp C_psn C_psp D_psn D_psp;

    % Electrolyte Current
    [F1_ien,F1_iep,F2_ien,F2_iep,F3_ien,F3_iep] = i_e_mats(p);
    p.F1_ien = F1_ien;
    p.F1_iep = F1_iep;
    p.F2_ien = F2_ien;
    p.F2_iep = F2_iep;
    p.F3_ien = F3_ien;
    p.F3_iep = F3_iep;

    clear F1_ien F1_iep F2_ien F2_iep F3_ien F3_iep;

    % Electrolyte Potential
    p.M1_pen_skel = sparse(diag(ones(p.Nxn-2,1),1) + diag(-ones(p.Nxn-2,1),-1));
    p.M1_pes_skel = sparse(diag(ones(p.Nxs-2,1),1) + diag(-ones(p.Nxs-2,1),-1));
    p.M1_pep_skel = sparse(diag(ones(p.Nxp-2,1),1) + diag(-ones(p.Nxp-2,1),-1));

    [M1_pe,M2_pe,M3_pe,M4_pe,C_pe] = phi_e_mats(p);
    p.M1_pe = M1_pe;
    p.M2_pe = M2_pe;
    p.M3_pe = M3_pe;
    p.M4_pe = M4_pe;
    p.C_pe = C_pe;

    clear M1_pe M2_pe M3_pe M4_pe C_pe

    %% CasADi Variable

    % Declare model variables
    x1 = SX.sym('x1', size(c_s_n,1)); % C_s_n
    x2 = SX.sym('x2', size(c_s_p,1)); % C_s_p
    x3 = SX.sym('x3', size(c_e,1));   % C_e
    x4 = SX.sym('x4', size(T10));     % T1
    x5 = SX.sym('x5', size(T20));     % T2

    x = [x1; x2; x3; x4; x5]; % State vector

    z1 = SX.sym('z1', size(phi_s_n,1)); % phi_s_n
    z2 = SX.sym('z2', size(phi_s_p,1)); % phi_s_p
    z3 = SX.sym('z3', size(i_en,1));    % i_en
    z4 = SX.sym('z4', size(i_ep,1));    % i_ep
    z5 = SX.sym('z5', size(phi_e,1));   % phi_e
    z6 = SX.sym('z6', size(jn,1));      % jn
    z7 = SX.sym('z7', size(jp,1));      % jp

    z = [z1; z2; z3; z4; z5; z6; z7]; % Algebraic variables

    % Input
    u = SX.sym('u'); 

    % Parameters
    % park=> parameter k

    %Sensitivity Params
    park = SX.zeros(length(sel_k),1); 
    for i=1:length(sel_k)
        park_tmp=cell_var{sel_k(i)};
        park(i)=park_tmp;
    end

    % Assign nominal sensitivity values
    park0 = zeros(length(sel_k),1);
    for i=1:length(sel_k)
        park0_tmp = SelecParam(i); % assign nominal values
        park0(i,1) = park0_tmp;
    end

    % Note, equilibrium params should be defined separately, i.e., eq_park
    % Total_park is sum of park & eq_park.
    % Not consider equilibrium sensitivity in init_cs.m.

    Total_park = park;
    Total_park0 = park0;

    Npark = size(Total_park,1);

    %% DAE Builder
%     [x_dot, g_, L ] = dae_dfn(x,z,u,p); % [SHP Change]
    [x_dot, g_, L, x_outs, z_outs, info_outs, param_outs ] = dae_dfn(x,z,u,p);

    %% Sensitivity Calculation (if turned on)

    % Initial condition
    x0_call = Function('x0_call',{park},{x0_nom},{'park'},{'result'});
    x0_init = full(x0_call(park0));
    x0 = full(x0_call(park0));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1. Calculate Jacobian automatically
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Jac_V0_param = jacobian(L,Total_park);
    V0_call = Function('V0_call',{park,x,u},{Jac_V0_param},{'park','x','u'},{'result'});
    S3_0_V0_params = full(V0_call(park0,x0_init,0));

    for j=1:length(sel_k)
        if sel_k(j) == 21 % c_e0 index
            Jac_x0_ce0 = jacobian(x0,p.c_e0);
            Jac_x0_ce0_call = Function('Jac_x0_ce0_call',{Total_park},{Jac_x0_ce0},{'Total_park'},{'result'});
            S1_0_x0_ce0 = full(Jac_x0_ce0_call(Total_park0));
        end
    end

    % Calculate Jacobian matrix (\frac{\partial f}{\partial x}}), A11
    f_x_jac = jacobian(x_dot,x);

    % Calculate Jacobian matrix (\frac{\partial f}{\partial z}}), A12
    f_z_jac = jacobian(x_dot,z);

    % Calculate Jacobian matrix (\frac{\partial f}{\partial \theta}}), B1
    f_theta_jac = jacobian(x_dot,Total_park);

    % Calculate Jacobian matrix (\frac{\partial g}{\partial x}}), A21
    g_x_jac = jacobian(g_,x);

    % Calculate Jacobian matrix (\frac{\partial g}{\partial z}}), A22
    g_z_jac = jacobian(g_,z);

    % Calculate Jacobian matrix (\frac{\partial g}{\partial \theta}}), B2
    g_theta_jac = jacobian(g_,Total_park);

    % Calculate Jacobian matrix (\frac{\partial h}{\partial x}}), C
    h_x_jac = jacobian(L,x);

    % Calculate Jacobian matrix (\frac{\partial h}{\partial z}}), D
    h_z_jac = jacobian(L,z);

    % Calculate Jacobian matrix (\frac{\partial h}{\partial \theta}}), E
    h_theta_jac = jacobian(L,Total_park);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 2. Build sensitivity equation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    S1_0 = zeros(size(x,1),Npark);
    for j=1:length(sel_k)
        pos = j;
        if sel_k(j) == 21 % c_e0 index
            S1_0(:,pos) = S1_0_x0_ce0;
        end
    end
    S2_0 = zeros(size(z,1),Npark);
    S3_0 = S3_0_V0_params;

    % Symbolic Sensitivity

    %%% S1_dot
    S1_blk = SX.sym('S1_blk',size(x,1)*Npark,1);
    S1 = reshape(S1_blk,size(x,1),Npark);

    S2_blk = SX.sym('S2_blk',size(z,1)*Npark,1);
    S2 = reshape(S2_blk,size(z,1),Npark);

    tmp = repmat({f_x_jac},Npark,1);
    f_x_blk = blkdiag(tmp{:});

    tmp = repmat({f_z_jac},Npark,1);
    f_z_blk = blkdiag(tmp{:});

    f_theta_blk = reshape(f_theta_jac,size(x,1)*Npark,1);

    S1_dot = f_x_blk*S1_blk + f_z_blk*S2_blk + f_theta_blk;

    %%% S2_
    tmp = repmat({g_x_jac},Npark,1);
    g_x_blk = blkdiag(tmp{:});

    tmp = repmat({g_z_jac},Npark,1);
    g_z_blk = blkdiag(tmp{:});

    g_theta_blk = reshape(g_theta_jac,size(z,1)*Npark,1);

    S2_ = g_x_blk*S1_blk + g_z_blk*S2_blk + g_theta_blk;

    %%% S3
    S3 = h_x_jac * S1 + h_z_jac * S2 + h_theta_jac;
    S3 = S3'; % change row-vector to column vector. 

    % IC for S1_blk, S2_blk
    S1_0_blk(:,1) = reshape(S1_0,size(x,1)*Npark,1);
    S2_0_blk(:,1) = reshape(S2_0,size(z,1)*Npark,1);

    %%% Explicitly define sensitivity values needed for simulation
    S1_sim_blk(:,1) = S1_0_blk(:,1);
    S2_sim_blk(:,1) = S2_0_blk(:,1);
    S3_sim(:,1) = S3_0';
    
    %% Integrator

    % Note: Sensitivity wouldnt converge/work for higher tolerance (10^-6) but
    % the DFN did, so separate CasADi integrator functions made for the DFN (F)
    % and Sensitivity DAEs (SS)
    del_t = p.delta_t;
    opts1 = struct('tf',del_t, 'abstol',1e-6,'reltol',1e-6);
    opts2 = struct('tf',del_t, 'abstol', 1e-1, 'reltol', 1e-1);

    dae = struct('x',x, 'z',z, 'p',[park;u], 'ode',x_dot, 'alg', g_, 'quad', L);
    F = integrator('F','idas', dae, opts1);

    if SensFlag == 1
        dae_s = struct('x',S1_blk, 'z', S2_blk, 'p',[park;x;z;u], 'ode', S1_dot, 'alg', S2_, 'quad', S3);
        SS = integrator('SS', 'idas', dae_s, opts2);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Build function for simulated states % [SHP Change]
    f_out = Function('f_out',{park,x,z,u},{x_outs},{'park','x','z','u'},{'f_out'});
    g_out = Function('g_out',{park,x,z,u},{z_outs},{'park','x','z','u'},{'g_out'});
    alg_out = Function('alg_out',{park,x,z,u},{info_outs},{'park','x','z','u'},{'alg_out'});
    par_out = Function('par_out',{park,x,z,u},{param_outs},{'park','x','z','u'},{'par_out'});
    f0 = full(f_out(park0,x0,z0,0));
    g0 = full(g_out(park0,x0,z0,0));
    a0 = full(alg_out(park0,x0,z0,0));
    p0 = full(par_out(park0,x0,z0,0));

    %% Indexing

    % index for x
    out_csn_idx = 1:(p.PadeOrder * (p.Nxn-1)); 
    out_csp_idx = (p.PadeOrder * (p.Nxn-1) + 1) : (p.PadeOrder * (p.Nxn-1) + p.PadeOrder * (p.Nxp-1));
    out_ce_idx = (p.PadeOrder * (p.Nxn-1) + p.PadeOrder * (p.Nxp-1) + 1) : (p.PadeOrder * (p.Nxn-1) + p.PadeOrder * (p.Nxp-1) + p.Nxn-1+p.Nxs-1+p.Nxp-1);
    out_T_idx = (p.PadeOrder * (p.Nxn-1) + p.PadeOrder * (p.Nxp-1) + p.Nxn-1+p.Nxs-1+p.Nxp-1 + 1);

%     disp('-------------')
%     disp('Check index in states X')
%     fprintf('size f0 = %d \n',size(f0,1))
%     fprintf(['cssn index:' int2str(out_csn_idx(1)) '...' int2str(out_csn_idx(end)) '\n'])
%     fprintf(['cssp index:' int2str(out_csp_idx(1)) '...' int2str(out_csp_idx(end)) '\n'])
%     fprintf(['cex index:' int2str(out_ce_idx(1)) '...' int2str(out_ce_idx(end)) '\n'])
%     fprintf(['T index:' int2str(out_T_idx(1)) '...' int2str(out_T_idx(end)) '\n'])

    % index for z

    out_phisn_idx = 1:p.Nxn-1; % 1:Nn
    out_phisp_idx = (p.Nxn-1+1) : (p.Nxn-1 + p.Nxp-1); % Nn+1 : Nnp
    out_ien_idx = (p.Nxn-1 + p.Nxp-1 +1) : (p.Nxn-1 + p.Nxp-1 + p.Nxn-1);
    out_iep_idx = (p.Nxn-1 + p.Nxp-1 + p.Nxn-1 +1) : (p.Nxn-1 + p.Nxp-1 + p.Nxn-1 + p.Nxp-1);
    out_phie_idx = (p.Nxn-1 + p.Nxp-1 + p.Nxn-1 + p.Nxp-1 + 1) : (p.Nxn-1 + p.Nxp-1 + p.Nxn-1 + p.Nxp-1 + (p.Nxn-1 + p.Nxs-1 + p.Nxp-1)+2); % Why 2? not 4?? used to plus 4 to consider boundary conditions.
    out_jn_idx = (p.Nxn-1 + p.Nxp-1 + p.Nxn-1 + p.Nxp-1 + (p.Nxn-1 + p.Nxs-1 + p.Nxp-1 + 2) + 1) : (p.Nxn-1 + p.Nxp-1 + p.Nxn-1 + p.Nxp-1 + (p.Nxn-1 + p.Nxs-1 + p.Nxp-1 + 2) + (p.Nxn-1));
    out_jp_idx = (p.Nxn-1 + p.Nxp-1 + p.Nxn-1 + p.Nxp-1 + (p.Nxn-1 + p.Nxs-1 + p.Nxp-1 + 2) + (p.Nxn-1) +1) : (p.Nxn-1 + p.Nxp-1 + p.Nxn-1 + p.Nxp-1 + (p.Nxn-1 + p.Nxs-1 + p.Nxp-1+ 2) + (p.Nxn-1) + (p.Nxp-1));
%     disp('-------------')
%     disp('Check index in states Z')
%     fprintf('size g0 = %d \n',size(g0,1))
%     fprintf(['phi_s_n index:' int2str(out_phisn_idx(1)) '...' int2str(out_phisn_idx(end)) '\n'])
%     fprintf(['phi_s_p index:' int2str(out_phisp_idx(1)) '...' int2str(out_phisp_idx(end)) '\n'])
%     fprintf(['i_en index:' int2str(out_ien_idx(1)) '...' int2str(out_ien_idx(end)) '\n'])
%     fprintf(['i_ep index:' int2str(out_iep_idx(1)) '...' int2str(out_iep_idx(end)) '\n'])
%     fprintf(['phi_e index:' int2str(out_phie_idx(1)) '...' int2str(out_phie_idx(end)) '\n'])
%     fprintf(['j_n index:' int2str(out_jn_idx(1)) '...' int2str(out_jn_idx(end)) '\n'])
%     fprintf(['j_p index:' int2str(out_jp_idx(1)) '...' int2str(out_jp_idx(end)) '\n'])

    % index for information
    out_cssn_idx = 1:p.Nxn-1 ;
    out_cssp_idx = (p.Nxn-1+1) : (p.Nxn-1 + p.Nxp-1);
    out_cex_idx = (p.Nxn-1 + p.Nxp-1 +1) : (p.Nxn-1 + p.Nxp-1 + (p.Nxn-1+p.Nxs-1+p.Nxp-1 +4));
    out_theta_avgn_idx = (p.Nxn-1 + p.Nxp-1 + (p.Nxn-1+p.Nxs-1+p.Nxp-1 +4) +1) : (p.Nxn-1 + p.Nxp-1 + (p.Nxn-1+p.Nxs-1+p.Nxp-1 +4) +(p.Nxn-1));
    out_theta_avgp_idx = (p.Nxn-1 + p.Nxp-1 + (p.Nxn-1+p.Nxs-1+p.Nxp-1 +4) +(p.Nxn-1) +1) : (p.Nxn-1 + p.Nxp-1 + (p.Nxn-1+p.Nxs-1+p.Nxp-1 +4) +(p.Nxn-1) +(p.Nxp-1));
    out_etan_idx = (p.Nxn-1 + p.Nxp-1 + (p.Nxn-1+p.Nxs-1+p.Nxp-1 +4) +(p.Nxn-1) +(p.Nxp-1) +1) : (p.Nxn-1 + p.Nxp-1 + (p.Nxn-1+p.Nxs-1+p.Nxp-1 +4) +(p.Nxn-1) +(p.Nxp-1) + (p.Nxn-1));
    out_etap_idx = (p.Nxn-1 + p.Nxp-1 + (p.Nxn-1+p.Nxs-1+p.Nxp-1 +4) +(p.Nxn-1) +(p.Nxp-1) + (p.Nxn-1) +1) : (p.Nxn-1 + p.Nxp-1 + (p.Nxn-1+p.Nxs-1+p.Nxp-1 +4) + (p.Nxn-1) +(p.Nxp-1) + (p.Nxn-1) + (p.Nxp-1));
    out_ce0n_idx = (p.Nxn-1 + p.Nxp-1 + (p.Nxn-1+p.Nxs-1+p.Nxp-1 +4) +(p.Nxn-1) +(p.Nxp-1) + (p.Nxn-1) + (p.Nxp-1) +1) : (p.Nxn-1 + p.Nxp-1 + (p.Nxn-1+p.Nxs-1+p.Nxp-1 +4) +(p.Nxn-1) +(p.Nxp-1) + (p.Nxn-1) + (p.Nxp-1) +1);
    out_ce0p_idx = (p.Nxn-1 + p.Nxp-1 + (p.Nxn-1+p.Nxs-1+p.Nxp-1 +4) +(p.Nxn-1) +(p.Nxp-1) + (p.Nxn-1) + (p.Nxp-1) +2) : (p.Nxn-1 + p.Nxp-1 + (p.Nxn-1+p.Nxs-1+p.Nxp-1 +4) +(p.Nxn-1) +(p.Nxp-1) + (p.Nxn-1) + (p.Nxp-1) +2);
    out_etasLn_idx = (p.Nxn-1 + p.Nxp-1 + (p.Nxn-1+p.Nxs-1+p.Nxp-1 +4) +(p.Nxn-1) +(p.Nxp-1) + (p.Nxn-1) + (p.Nxp-1) +3) : (p.Nxn-1 + p.Nxp-1 + (p.Nxn-1+p.Nxs-1+p.Nxp-1 +4) +(p.Nxn-1) +(p.Nxp-1) + (p.Nxn-1) + (p.Nxp-1) +3);
    out_Volt_idx = (p.Nxn-1 + p.Nxp-1 + (p.Nxn-1+p.Nxs-1+p.Nxp-1 +4) +(p.Nxn-1) +(p.Nxp-1) + (p.Nxn-1) + (p.Nxp-1) +4) : (p.Nxn-1 + p.Nxp-1 + (p.Nxn-1+p.Nxs-1+p.Nxp-1 +4) +(p.Nxn-1) +(p.Nxp-1) + (p.Nxn-1) + (p.Nxp-1) +4);
    out_nLis_idx = (p.Nxn-1 + p.Nxp-1 + (p.Nxn-1+p.Nxs-1+p.Nxp-1 +4) +(p.Nxn-1) +(p.Nxp-1) + (p.Nxn-1) + (p.Nxp-1) +5) : (p.Nxn-1 + p.Nxp-1 + (p.Nxn-1+p.Nxs-1+p.Nxp-1 +4) +(p.Nxn-1) +(p.Nxp-1) + (p.Nxn-1) + (p.Nxp-1) +5);
    out_nLie_idx = (p.Nxn-1 + p.Nxp-1 + (p.Nxn-1+p.Nxs-1+p.Nxp-1 +4) +(p.Nxn-1) +(p.Nxp-1) + (p.Nxn-1) + (p.Nxp-1) +6) : (p.Nxn-1 + p.Nxp-1 + (p.Nxn-1+p.Nxs-1+p.Nxp-1 +4) +(p.Nxn-1) +(p.Nxp-1) + (p.Nxn-1) + (p.Nxp-1) +6);

%     disp('-------------')
%     disp('Check index in simulation info. (algebraics)')
%     fprintf('size a0 = %d \n',size(a0,1))
%     fprintf(['cssn index:' int2str(out_cssn_idx(1)) '...' int2str(out_cssn_idx(end)) '\n'])
%     fprintf(['cssp index:' int2str(out_cssp_idx(1)) '...' int2str(out_cssp_idx(end)) '\n'])
%     fprintf(['cex index:' int2str(out_cex_idx(1)) '...' int2str(out_cex_idx(end)) '\n'])
%     fprintf(['theta_avgn index:' int2str(out_theta_avgn_idx(1)) '...' int2str(out_theta_avgn_idx(end)) '\n'])
%     fprintf(['theta_avgp index:' int2str(out_theta_avgp_idx(1)) '...' int2str(out_theta_avgp_idx(end)) '\n'])
%     fprintf(['etan index:' int2str(out_etan_idx(1)) '...' int2str(out_etan_idx(end)) '\n'])
%     fprintf(['etap index:' int2str(out_etap_idx(1)) '...' int2str(out_etap_idx(end)) '\n'])
%     fprintf(['ce0n index:' int2str(out_ce0n_idx(1)) '...' int2str(out_ce0n_idx(end)) '\n'])
%     fprintf(['ce0p index:' int2str(out_ce0p_idx(1)) '...' int2str(out_ce0p_idx(end)) '\n'])
%     fprintf(['etasLn index:' int2str(out_etasLn_idx(1)) '...' int2str(out_etasLn_idx(end)) '\n'])
%     fprintf(['Volt index:' int2str(out_Volt_idx(1)) '...' int2str(out_Volt_idx(end)) '\n'])
%     fprintf(['nLis index:' int2str(out_nLis_idx(1)) '...' int2str(out_nLis_idx(end)) '\n'])
%     fprintf(['nLie index:' int2str(out_nLie_idx(1)) '...' int2str(out_nLie_idx(end)) '\n'])

    % index for state-dep param
    out_Den_idx = 1:(p.Nxn-1) ;
    out_dDen_idx = (p.Nxn-1 +1) : (p.Nxn-1 + p.Nxn-1);
    out_Des_idx = (p.Nxn-1 + p.Nxn-1 +1) : (p.Nxn-1 + p.Nxn-1 + p.Nxs-1);
    out_dDes_idx = (p.Nxn-1 + p.Nxn-1 + p.Nxs-1 +1) : (p.Nxn-1 + p.Nxn-1 + p.Nxs-1 + p.Nxs-1);
    out_Dep_idx = (p.Nxn-1 + p.Nxn-1 + p.Nxs-1 + p.Nxs-1 +1) : (p.Nxn-1 + p.Nxn-1 + p.Nxs-1 + p.Nxs-1 + p.Nxp -1);
    out_dDep_idx = (p.Nxn-1 + p.Nxn-1 + p.Nxs-1 + p.Nxs-1 + p.Nxp-1 +1) : (p.Nxn-1 + p.Nxn-1 + p.Nxs-1 + p.Nxs-1 + p.Nxp-1 + p.Nxp-1);

%     disp('Check index in simulation info. (params)')
%     fprintf('size p0 = %d \n',size(p0,1))
%     fprintf(['Den0 index:' int2str(out_Den_idx(1)) '...' int2str(out_Den_idx(end)) '\n']);
%     fprintf(['dDen0 index:' int2str(out_dDen_idx(1)) '...' int2str(out_dDen_idx(end)) '\n']);
%     fprintf(['Des0 index:' int2str(out_Des_idx(1)) '...' int2str(out_Des_idx(end)) '\n']);
%     fprintf(['dDes0 index:' int2str(out_dDes_idx(1)) '...' int2str(out_dDes_idx(end)) '\n']);
%     fprintf(['Dep0 index:' int2str(out_Dep_idx(1)) '...' int2str(out_Dep_idx(end)) '\n']);
%     fprintf(['dDep0 index:' int2str(out_dDep_idx(1)) '...' int2str(out_dDep_idx(end)) '\n']);

    %% Simulation

    % pre-allocate
    v_sim = zeros(NT,1);
    T1_sim = zeros(NT,1);
    T2_sim = zeros(NT,1);
    sens = zeros(NT,np);
    S3_sim = zeros(np,NT);
    
    % Initialize
    x_sim(:,1) = x0;
    z_sim(:,1) = z0;
    v_sim(1) = V0;

    % [SHP Change]
    % save x
    csn_sim(:,1) = f0(out_csn_idx);
    csp_sim(:,1) = f0(out_csp_idx);
    ce_sim(:,1) = f0(out_ce_idx);
    cen_sim(:,1) = ce_sim(1:(p.Nxn-1),1);
    ces_sim(:,1) = ce_sim((p.Nxn-1 +1):(p.Nxn-1 + p.Nxs-1),1);
    cep_sim(:,1) = ce_sim((p.Nxn-1 + p.Nxs-1 +1):(p.Nxn-1 + p.Nxs-1 + p.Nxp-1),1);
    T1_sim(1) = f0(out_T_idx);
    T2_sim(1) = f0(end);
    
    % save z
    phi_s_n_sim(:,1) = g0(out_phisn_idx);
    phi_s_p_sim(:,1) = g0(out_phisp_idx);
    ien_sim(:,1) = g0(out_ien_idx);
    iep_sim(:,1) = g0(out_iep_idx);
    phie_sim(:,1) = g0(out_phie_idx);
    jn_sim(:,1) = g0(out_jn_idx);
    jp_sim(:,1) = g0(out_jp_idx);

    % save alg. states
    cssn_sim(:,1) = a0(out_cssn_idx);
    cssp_sim(:,1) = a0(out_cssp_idx);
    cex_sim(:,1) = a0(out_cex_idx);
    theta_avgn_sim(:,1) = a0(out_theta_avgn_idx);
    theta_avgp_sim(:,1) = a0(out_theta_avgp_idx);
    etan_sim(:,1) = a0(out_etan_idx);
    etap_sim(:,1) = a0(out_etap_idx);
    ce0n_sim(:,1) = a0(out_ce0n_idx);
    ce0p_sim(:,1) = a0(out_ce0p_idx);
    etasLn_sim(:,1) = a0(out_etasLn_idx);
    volt_sim(:,1) = a0(out_Volt_idx);
    nLis_sim(:,1) = a0(out_nLis_idx);
    nLie_sim(:,1) = a0(out_nLie_idx);

    % save state-dependent param
    Den0_sim(:,1) = p0(out_Den_idx);
    dDen0_sim(:,1) = p0(out_dDen_idx);
    Des0_sim(:,1) = p0(out_Des_idx);
    dDes0_sim(:,1) = p0(out_dDes_idx);
    Dep0_sim(:,1) = p0(out_Dep_idx);
    dDep0_sim(:,1) = p0(out_dDep_idx);

%     disp('---------------------------------')
%     disp('Running')
    
    % Simulate DFN & Sensitivity Equations
    try
        for k=1:(NT-1)
            if (mod(k,500) == 0)
                fprintf('Iter:%d, v_sim: %f\n',k,v_sim(k));
            end

            Cur = I(k+1);

            % DFN
            Fk = F('x0',x_sim(:,k),'z0',z_sim(:,k),'p',[park0;Cur]);
            x_sim(:,k+1) = full(Fk.xf);
            z_sim(:,k+1) = full(Fk.zf);
            v_sim(k+1) = full(Fk.qf)/p.delta_t;

            f0 = full(f_out(park0,x_sim(:,k+1),z_sim(:,k+1),Cur));
            g0 = full(g_out(park0,x_sim(:,k+1),z_sim(:,k+1),Cur));
            a0 = full(alg_out(park0,x_sim(:,k+1),z_sim(:,k+1),Cur));
            p0 = full(par_out(park0,x_sim(:,k+1),z_sim(:,k+1),Cur));

            % save x
            csn_sim(:,k+1) = f0(out_csn_idx);
            csp_sim(:,k+1) = f0(out_csp_idx);
            ce_sim(:,k+1) = f0(out_ce_idx);
            cen_sim(:,k+1) = ce_sim(1:(p.Nxn-1),k+1);
            ces_sim(:,k+1) = ce_sim((p.Nxn-1 +1):(p.Nxn-1 + p.Nxs-1),k+1);
            cep_sim(:,k+1) = ce_sim((p.Nxn-1 + p.Nxs-1 +1):(p.Nxn-1 + p.Nxs-1 + p.Nxp-1),k+1);

            T1_sim(k+1) = f0(out_T_idx);
            T2_sim(k+1) = f0(end);

            % save z
            phi_s_n_sim(:,k+1) = g0(out_phisn_idx);
            phi_s_p_sim(:,k+1) = g0(out_phisp_idx);
            ien_sim(:,k+1) = g0(out_ien_idx);
            iep_sim(:,k+1) = g0(out_iep_idx);
            phie_sim(:,k+1) = g0(out_phie_idx);
            jn_sim(:,k+1) = g0(out_jn_idx);
            jp_sim(:,k+1) = g0(out_jp_idx);

            % save alg.
            cssn_sim(:,k+1) = a0(out_cssn_idx);
            cssp_sim(:,k+1) = a0(out_cssp_idx);
            cex_sim(:,k+1) = a0(out_cex_idx);
            theta_avgn_sim(:,k+1) = a0(out_theta_avgn_idx);
            theta_avgp_sim(:,k+1) = a0(out_theta_avgp_idx);
            etan_sim(:,k+1) = a0(out_etan_idx);
            etap_sim(:,k+1) = a0(out_etap_idx);
            ce0n_sim(:,k+1) = a0(out_ce0n_idx);
            ce0p_sim(:,k+1) = a0(out_ce0p_idx);
            etasLn_sim(:,k+1) = a0(out_etasLn_idx);
            volt_sim(:,k+1) = a0(out_Volt_idx);
            nLis_sim(:,k+1) = a0(out_nLis_idx);
            nLie_sim(:,k+1) = a0(out_nLie_idx);

            % save param.
            Den0_sim(:,k+1) = p0(out_Den_idx);
            dDen0_sim(:,k+1) = p0(out_dDen_idx);
            Des0_sim(:,k+1) = p0(out_Des_idx);
            dDes0_sim(:,k+1) = p0(out_dDes_idx);
            Dep0_sim(:,k+1) = p0(out_Dep_idx);
            dDep0_sim(:,k+1) = p0(out_dDep_idx);

            % Step SOC forward
            SOC(k+1) = (mean(c_avg_n(:,k+1)) - cn_low) / (cn_high - cn_low);

            if SensFlag == 1 % Calculate sensitivity flag == true
                % Sensitivity eqns
                Sk = SS('x0',S1_sim_blk(:,k),'z0',S2_sim_blk(:,k),'p',[park0;x_sim(:,k);z_sim(:,k);Cur]);
                S1_sim_blk(:,k+1) = full(Sk.xf);
                S2_sim_blk(:,k+1) = full(Sk.zf);
                S3_sim(:,k+1) = full(Sk.qf)/p.delta_t;
            end

             % Check voltage constraints and exit if violated
            if v_sim(k+1) <= p.volt_min
                fprintf('Min voltage is reached at %d iteration. Stopping simulation and concatenating data. \n',k);
                v_sim = concatenate_data(v_sim);
                %%%%%%%%%%%%%    ZTG Note: should the rest of the gradient be zero or the last gradient value? 
                sens(1:length(S3_sim),:) = S3_sim'; 
                varargout{1} = sens;
                alg_states = [];
                return
            end

            if v_sim(k+1) >= p.volt_max
                fprintf('Max voltage is reached at %d iteration. Stopping simulation and concatenating data. \n',k);
                v_sim = concatenate_data(v_sim);
                %%%%%%%%%%%%%    ZTG Note: should the rest of the gradient be zero or the last gradient value? 
                sens(1:length(S3_sim),:) = S3_sim';
                varargout{1} = sens;
                alg_states = [];
                return 
            end
        end
    catch e
        errorMessage = sprintf('%s',getReport( e, 'extended', 'hyperlinks', 'on' ))
        fprintf('CasADi error. Stopping simulation and concatenating data. \n');
        v_sim = concatenate_data(v_sim);
        %%%%%%%%%%%%%    ZTG Note: should the rest of the gradient be zero or the last gradient value? 
        sens(1:length(S3_sim),:) = S3_sim';
        varargout{1} = sens;

        alg_states.T1_sim = T1_sim;
        alg_states.T2_sim = T2_sim;
        alg_states.phi_s_n_sim = phi_s_n_sim;
        alg_states.phi_s_p_sim = phi_s_p_sim;
        alg_states.phie_sim = phie_sim;
        alg_states.ien_sim = ien_sim;
        alg_states.iep_sim = iep_sim;
        alg_states.jn_sim = jn_sim;
        alg_states.jp_sim = jp_sim;
        alg_states.cssn_sim = cssn_sim;
        alg_states.cssp_sim = cssp_sim;
        alg_states.cex_sim = cex_sim;
        alg_states.theta_avgn_sim = theta_avgn_sim;
        alg_states.theta_avgp_sim = theta_avgp_sim;
        alg_states.etan_sim = etan_sim;
        alg_states.etap_sim = etap_sim;
        
        save('casadi_debug.mat','v_sim','alg_states','Current_exp', 'Time_exp')
        return
    end

    %% Collect Outputs
%     alg_states.csn_sim = csn_sim;
%     alg_states.csp_sim = csp_sim;
%     alg_states.ce_sim = ce_sim;
    alg_states.T1_sim = T1_sim;
    alg_states.T2_sim = T2_sim;
%     alg_states.phi_s_n_sim = phi_s_n_sim;
%     alg_states.phi_s_p_sim = phi_s_p_sim;
%     alg_states.ien_sim = ien_sim;
%     alg_states.iep_sim = iep_sim;
%     alg_states.phie_sim = phie_sim;
%     alg_states.jn_sim = jn_sim;
%     alg_states.jp_sim = jp_sim;
    alg_states.cssn_sim = cssn_sim;
    alg_states.cssp_sim = cssp_sim;
%     alg_states.cex_sim = cex_sim;
%     alg_states.theta_avgn_sim = theta_avgn_sim;
%     alg_states.theta_avgp_sim = theta_avgp_sim;
    alg_states.etan_sim = etan_sim;
    alg_states.etap_sim = etap_sim;
%     alg_states.ce0n_sim = ce0n_sim;
%     alg_states.ce0p_sim = ce0p_sim;
%     alg_states.etasLn_sim = etasLn_sim;
%     alg_states.volt_sim = volt_sim;
%     alg_states.nLis_sim = nLis_sim;
%     alg_states.nLie_sim = nLie_sim;
%     alg_states.SOC = SOC;
%     alg_states.p = p;
%     alg_states.Den0_sim = Den0_sim;
%     alg_states.dDen0_sim = dDen0_sim;
%     alg_states.Des0_sim = Des0_sim;
%     alg_states.dDes0_sim = dDes0_sim;
%     alg_states.Dep0_sim = Dep0_sim;
%     alg_states.dDep0_sim = dDep0_sim;
%     alg_states.cen_sim = cen_sim;
%     alg_states.ces_sim = ces_sim;
%     alg_states.cep_sim = cep_sim;

    %% Return
    if SensFlag == 1
        sens = S3_sim'; % Nt-by-Np matrix
        varargout{1} = sens;
    end
end

