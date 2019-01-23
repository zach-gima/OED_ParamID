%% "Nominal" Parameters from Literature Used to Initialize ParamID by Saehong Park in Bosch-UC Berkeley Yr. 1
% Corresponds to NCA cell (param/params_NCA_final.m)
% By: Zach Gima 2018-2-22
%% May/may not need to add in the 2-state thermal dynamics parameters into this

% Nominal_param = [3.90e-14;    %1.D_s_n
%                  1.00e-13;    %2.D_s_p
%                  10.9e-06;   %3.R_s_n
%                  10.9e-06;   %4.R_s_p
%                  nan;        %5.eps_s_n (x), Equil. Struct.
%                  nan;        %6.eps_s_p (x), Equil. Struct.
%                  100;        %7.sig_n
%                  100;        %8.sig_p
%                  1;          %9.D_e
%                  0.3;        %10.eps_e_n
%                  0.5;        %11.eps_e_s
%                  0.3;        %12.eps_e_p
%                  1;          %13.Kappa
%                  0.363;       %14.t_plus
%                  1;          %15.d_activity
%                  7.5e-04;    %16.k_n0
%                  2.3e-03;    %17.k_p0
%                  5e-04;      %18.R_f_n
%                  1e-03;      %19.R_f_p
%                  nan;        %20.n_Li_s (x), Equil. Struct.
%                  1e3;         %21.ce0
%                  36.63e3;       %22.E.Dsn
%                  47.98e3;       %23.E.Dsp
%                  53.4e3;       %24.E.kn
%                  39.57e3        %25.E.kp                 
%                  ]; 
             
%%% Experimentally Identified Parameters 2019-1-15
expID_param =   [1.05000000000000e-12;    %1.D_s_n
                 8.12920311552651e-13;    %2.D_s_p
                 9.72819296281088e-05;   %3.R_s_n
                 5.32697934211608e-05;   %4.R_s_p
                 nan;        %5.eps_s_n (x), Equil. Struct.
                 nan;        %6.eps_s_p (x), Equil. Struct.
                 100;        %7.sig_n
                 100;        %8.sig_p
                 1.50000000000000;          %9.D_e
                 0.3;        %10.eps_e_n
                 0.500000000000000;        %11.eps_e_s
                 0.272256850211865;        %12.eps_e_p
                 1;          %13.Kappa
                 0.38;       %14.t_plus 0.36 originally; 0.38 V-R Relationship
                 1.50000000000000;          %15.d_activity
                 7.5e-04;    %16.k_n0
                 2.3e-03;    %17.k_p0
                 0.000589404729676806;      %18.R_f_n
                 1e-03;      %19.R_f_p
                 nan;        %20.n_Li_s (x), Equil. Struct.
                 500;         %21.ce0
                 41417.9392555510;       %22.E.Dsn
                 28459.7916534643;       %23.E.Dsp
                 67500;       %24.E.kn
                 39.57e3        %25.E.kp                 
                 ];              