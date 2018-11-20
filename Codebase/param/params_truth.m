%% "Truth" Parameters Identified by Saehong Park in Bosch-UC Berkeley Yr. 1 (as of 2018-2-22)
% Corresponds to NCA cell (param/params_NCA_final.m)
% By: Zach Gima 2018-2-22
 
%% May/may not need to add in the 2-state thermal dynamics parameters into this

%%% LIT REVIEW THERMAL PARAMS
% truth_param = [2.634e-14;    %1.D_s_n
%                  6.625e-14;    %2.D_s_p
%                  2.0235e-05;   %3.R_s_n
%                  1.7163e-05;   %4.R_s_p
%                  0.5438;%nan;        %5.eps_s_n (x), Equil. Struct.
%                  0.6663;%nan;        %6.eps_s_p (x), Equil. Struct.
%                  500;        %7.sig_n
%                  500;        %8.sig_p
%                  1.195;          %9.D_e
%                  0.289;        %10.eps_e_n
%                  0.468;        %11.eps_e_s
%                  0.307;        %12.eps_e_p
%                  1.398;          %13.Kappa
%                  0.36;       %14.t_plus
%                  0.573;          %15.d_activity
%                  7.5e-05;    %16.k_n0
%                  2.3e-04;    %17.k_p0
%                  8.719e-05;      %18.R_f_n
%                  4.619e-04;      %19.R_f_p
%                  0.1406;%nan;        %20.n_Li_s (x), Equil. Struct.
%                  1.5e3;         %21.ce0
%                  36.63e3;       %22.E.Dsn
%                  47.98e3;       %23.E.Dsp
%                  53.4e3;       %24.E.kn
%                  39.57e3        %25.E.kp
%                  ];  % assign nominal vlaues that selects sensitivity

%%% ZHANG THERMAL PARAMS
truth_param = [2.634e-14;    %1.D_s_n
                 6.625e-14;    %2.D_s_p
                 2.0235e-05;   %3.R_s_n
                 1.7163e-05;   %4.R_s_p
                 0.5438;%nan;        %5.eps_s_n (x), Equil. Struct.
                 0.6663;%nan;        %6.eps_s_p (x), Equil. Struct.
                 500;        %7.sig_n
                 500;        %8.sig_p
                 1.195;          %9.D_e
                 0.289;        %10.eps_e_n
                 0.468;        %11.eps_e_s
                 0.307;        %12.eps_e_p
                 1.398;          %13.Kappa
                 0.38;       %14.t_plus %0.36 originally; 0.38 V-R Relationship
                 0.573;          %15.d_activity
                 7.5e-05;    %16.k_n0
                 2.3e-04;    %17.k_p0
                 8.719e-05;      %18.R_f_n
                 4.619e-04;      %19.R_f_p
                 0.1406;%nan;        %20.n_Li_s (x), Equil. Struct.
                 1.5e3;         %21.ce0
                 36.63e3;       %22.E.Dsn % 2018-9-20 currently nominal values from lit. [ZTG]
                 47.98e3;       %23.E.Dsp % 2018-9-20 currently nominal values from lit. [ZTG]
                 53.4e3;       %24.E.kn % 2018-9-20 currently nominal values from lit. [ZTG]
                 39.57e3        %25.E.kp % 2018-9-20 currently nominal values from lit. [ZTG]
                 ];  % assign nominal vlaues that selects sensitivity             