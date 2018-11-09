%% Params for Electrochemical Model
%  Created Apr 26, 2016 by Saehong Park
%
%   Most parameters are from NCA-18650
%   D_e is from Capiglia et al, for c_e = 1000 mol/m^3

% Ref #1
% Uddin et al. (Batteries 2016) (charge degrad.)								
% https://www.mendeley.com/share/document/invite/638008b478/?utm_medium=email&utm_source=transactional&utm_campaign=share%2Finvitation-document								

% Ref #2
% T.R. Ashwin (Power Sources 2016) (mod. El. Par.)
% https://www.mendeley.com/share/document/invite/a3d3f3d0b7/?utm_medium=email&utm_source=transactional&utm_campaign=share%2Finvitation-document

% Ref #3
% Bernadi et al (JPS 2011)
% http://www.sciencedirect.com/science/article/pii/S0378775310011468

%% Geometric Params
% Thickness of each layer
p.L_n = 79e-6;     % Thickness of negative electrode [m]
p.L_s = 80e-6;     % Thickness of separator [m]
p.L_p = 61.5e-6;     % Thickness of positive electrode [m]

L_ccn = 25e-6;    % Thickness of negative current collector [m]
L_ccp = 25e-6;    % Thickness of negative current collector [m]

% Particle Radii
p.R_s_n = 2.0235345e-05;%(True)   % Radius of solid particles in negative electrode [m]
p.R_s_p = 1.7162987e-05;   % Radius of solid particles in positive electrode [m]

% Volume fractions
p.epsilon_s_n = 0.543889597565723;%0.6987;      % Volume fraction in solid for neg. electrode
p.epsilon_s_p = 0.666364981170368;%0.5291;      % Volume fraction in solid for pos. electrode

p.epsilon_e_n = 2.8908315e-01;   % Volume fraction in electrolyte for neg. electrode
p.epsilon_e_s = 4.6798515e-01;   % Volume fraction in electrolyte for separator
p.epsilon_e_p = 3.0676745e-01;   % Volume fraction in electrolyte for pos. electrode

% make element to caclulate phi_{s} by Saehong Park 
p.epsilon_f_n = 1 - p.epsilon_s_n - p.epsilon_e_n;  % Volume fraction of filler in neg. electrode
p.epsilon_f_p = 1 - p.epsilon_s_p - p.epsilon_e_p;  % Volume fraction of filler in pos. electrode

epsilon_f_n = p.epsilon_f_n;  % Volume fraction of filler in neg. electrode
epsilon_f_p = p.epsilon_f_p;  % Volume fraction of filler in pos. electrode


% Specific interfacial surface area
p.a_s_n = 3*p.epsilon_s_n / p.R_s_n;  % Negative electrode [m^2/m^3]
p.a_s_p = 3*p.epsilon_s_p / p.R_s_p;  % Positive electrode [m^2/m^3]

% % % % % % % % % % Mass densities
% % % % % % % % % rho_sn = 1800;    % Solid phase in negative electrode [kg/m^3]
% % % % % % % % % rho_sp = 5010;    % Solid phase in positive electrode [kg/m^3]
% % % % % % % % % rho_e =  1324;    % Electrolyte [kg/m^3]
% % % % % % % % % rho_f = 1800;     % Filler [kg/m^3]
% % % % % % % % % rho_ccn = 8954;   % Current collector in negative electrode
% % % % % % % % % rho_ccp = 2707;   % Current collector in positive electrode
% % % % % % % % % 
% % % % % % % % % % Compute cell mass [kg/m^2]
% % % % % % % % % m_n = p.L_n * (rho_e*p.epsilon_e_n + rho_sn*p.epsilon_s_n + rho_f*epsilon_f_n);
% % % % % % % % % m_s = p.L_s * (rho_e*p.epsilon_e_n);
% % % % % % % % % m_p = p.L_p * (rho_e*p.epsilon_e_p + rho_sp*p.epsilon_s_p + rho_f*epsilon_f_p);
% % % % % % % % % m_cc = rho_ccn*L_ccn + rho_ccp*L_ccp;

% Lumped density [kg/m^2]
% % % % % % % % % p.rho_avg = m_n + m_s + m_p + m_cc;

%% Transport Params
% Diffusion coefficient in solid
p.D_s_n0 = 2.6340210e-14;%3.9e-14(True)   % Diffusion coeff for solid in neg. electrode, [m^2/s]
p.D_s_p0 = 6.6259080e-14; % Diffusion coeff for solid in pos. electrode, [m^2/s]

% Diffusion coefficient for CVM 
p.D_s_n_Factor = 1;
p.D_s_p_Factor = 1;

% Electrolyte Perturbation Factor
p.ElecFactorK = 1.398236e+00; % used for ElectrolyteConductivity
p.ElecFactorD = 1.195133e+00; % used for Electrolyte Diffusivity
p.ElecFactorDA = 5.73455e-01;% used for dactivity

% Diffusional conductivity in electrolyte
p.dactivity = 0;

p.brug = 1.5;       % Bruggeman porosity

% Conductivity of solid
p.sig_n = 5.000000e+02;    % Conductivity of solid in neg. electrode, [1/Ohms*m]
p.sig_p = 5.000000e+02;%(True)    % Conductivity of solid in pos. electrode, [1/Ohms*m]

% Miscellaneous
p.t_plus = 3.600000e-01;       % Transference number
p.Faraday = 96487;    % Faraday's constant, [Coulumbs/mol]
p.Area = 0.0744;           % Electrode current collector area [m^2]

%% Kinetic Params
p.R = 8.314472;       % Gas constant, [J/mol-K]

p.alph = 0.5;         % Charge transfer coefficients

p.R_f_n = 8.7191076e-05;       % Resistivity of SEI layer, [Ohms*m^2]
p.R_f_p = 4.6199623e-04;       % Resistivity of SEI layer, [Ohms*m^2]
p.R_c = 2.5e-03;%5.1874e-05/p.Area; % Contact Resistance/Current Collector Resistance, [Ohms-m^2]

% Nominal Reaction rates
p.k_n0 = 7.500000e-05;  % Reaction rate in neg. electrode, [(A/m^2)*(mol^3/mol)^(1+alpha)]
p.k_p0 = 2.300000e-04; % Reaction rate in pos. electrode, [(A/m^2)*(mol^3/mol)^(1+alpha)]

%% Thermodynamic Params

% Thermal dynamics
p.C_p = 2000;   % Heat capacity, [J/kg-K]
p.h = 0.36;   % Heat transfer coefficient, [W/K-m^2] 0

% Ambient Temperature
p.T_amb = 298.15; % [K]

% Activation Energies
% Taken from Zhang et al (2014) [Harbin]
% http://dx.doi.org/10.1016/j.jpowsour.2014.07.110
% All units are [J/mol]
p.E.kn = 37.48e3;
p.E.kp = 39.57e3;
p.E.Dsn = 42.77e3;
p.E.Dsp = 18.55e3;
p.E.De = 37.04e3;
p.E.kappa_e = 34.70e3;

% Reference temperature
p.T_ref = 298.15; %[K]

% % % Heat transfer parameters
% % % Taken from Zhang et al (2014) [Harbin]
% % % http://dx.doi.org/10.1016/j.jpowsour.2014.07.110
% % p.C1 = 62.7;    % [J/K]
% % p.C2 = 4.5;     % [J/K]
% % p.h12 = 10; % [W/K]
% % p.h2a = 21.45;  % [W/K]

% Heat transfer parameters
% Taken from Zhang et al (2014) [Harbin]
% http://dx.doi.org/10.1016/j.jpowsour.2014.07.110
p.C1 = 62.7/p.Area;    % [J/K]
p.C2 = 4.5/p.Area;     % [J/K]
p.h12 = (1/10)/p.Area; % [W/K]
p.h2a = (1/21.45)/p.Area;  % [W/K]

% % % Heat transfer parameters
% % % Taken from Hector's disseratation
% % p.C1 = 62.7/p.Area;    % [J/(m^2-K)]
% % p.C2 = 4.5/p.Area;     % [J/(m^2-K)]
% % p.h12 = (1/1.94)/p.Area; % [W/(m^2-K)]
% % p.h2a = 1.5*(1/3.08)/p.Area;  % [W/(m^2-K)]



%% Concentrations
% Maxima based on DUALFOIL 
% line 588 in DUALFOIL Fortran code

p.c_s_n_max = 3.71e+04;   % Max concentration in anode, [mol/m^3]

%p.c_s_n_max = 3.6e3 * 372 * 1800 / p.Faraday;   % Max concentration in anode, [mol/m^3]


p.c_s_p_max = 5.10e+04;    % Max concentration in cathode, [mol/m^3]

%p.c_s_p_max = 3.6e3 * 274 * 5010 / p.Faraday;    % Max concentration in cathode, [mol/m^3]

p.n_Li_s = 0.1406; %2.5; %2.781;        % Total moles of lithium in solid phase [mol]
p.c_e0 = 1.500000e+03;              % Fixed electrolyte concentration for SPM, [mol/m^3]

%% Cutoff voltages
p.volt_max = 4.2; %4.1113; %4.7;
p.volt_min = 2.5; %2.6;

%% Discretization parameters

% Pade Order
p.PadeOrder = 3;

% Finite difference points along x-coordinate
p.Nxn = 70;
p.Nxs = 35;
p.Nxp = 70;
p.Nx = p.Nxn+p.Nxs+p.Nxp;

p.delta_x_n = 1 / p.Nxn;
p.delta_x_s = 1 / p.Nxs;
p.delta_x_p = 1 / p.Nxp;

