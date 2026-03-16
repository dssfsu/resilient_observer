%% Resilient Observer Initializer
%   - Loads bus data
%   - Defines network Laplacian
%   - Defines discretize system matrices
%   - Resilient Observer
%   - Attack Model
%   - Bad Data Detection
%   - Auxiliary Model

% Jan. 2020
% Olugbenga Moses Anubi
% Carlos A. Wong
% Satish Vedula

% March 2026
% Olugbenga Moses Anubi
% Koto Andrew Omiloli

%% Parameter Setup
clear variables
close all
clc

addpath('Resilient_Optimizer')
addpath('Attack_Generators')

%% Data Loading (MATPOWER)

%=== MATPOWER global indices ==============================================
% column indices for bus data
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;

% column indices for branch data
[F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, ...
    RATE_C, TAP, SHIFT, BR_STATUS, PF, QF, PT, QT, MU_SF, MU_ST, ...
    ANGMIN, ANGMAX, MU_ANGMIN, MU_ANGMAX] = idx_brch;

% column indices for gen data
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;

% ========= CREATING CASE STRUCTURE =============================

mpc = loadcase('case39');

baseMVA     = mpc.baseMVA;
bus         = mpc.bus;
gen         = mpc.gen;
branch      = mpc.branch;

%% Dimension Variables
n_gen       = size(gen,1); % number of generators
n_bus       = size(bus,1); % number of buses
n_states    = 2*n_gen;
p           = n_states + n_bus;
n_load      = n_bus-n_gen;
n_input     = n_gen + n_load;  % number of measurements


gen_bus_idx = mpc.gen(:,GEN_BUS);
load_bus_idx = setdiff(1:n_bus, gen_bus_idx);     % remaining load buses

% Slack bus
idx_slack = bus(bus(:,BUS_TYPE)==3,BUS_I);  % slack bus number.
slack_mask = gen(:,GEN_BUS)==idx_slack;

%% Network Laplacian
[Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);

% Step1: get conductance and subsceptance matrix
G = real(Ybus);   % Conductance matrix
B = imag(Ybus);   % Susceptance matrix

%Step 2: get voltage magnitude
Vm = bus(:,VM);     % Voltage magnitude

%Step 3: Get Laplacian (Linearizing about equilibrium points (Vm, Thetas))
Lap = diag(Vm) * ((diag(B * Vm) - B * diag(Vm)));

%Step 4: Extracting Laplacian sub matrices
L_gg = Lap(gen_bus_idx, gen_bus_idx);        % generator-generator block
L_gl = Lap(gen_bus_idx, load_bus_idx);       % generator-load
L_lg = L_gl.';                       % load-generator
L_ll = Lap(load_bus_idx, load_bus_idx);      % load-load

%Step 5: Ssanity check
disp('Sanity Check graph Laplacian: Positive semi-definite?')
disp(eig(L_gg-(L_gl*(L_ll\L_lg))).')

%% Inertia Matrix
% Reference data
f_ref = 60;             % reference frequency in Hz
omega_r = 2*pi*f_ref;   % reference frequency in rad/sec

% Generator inertia constant and damping
H_range = [3 8];  % range of inertia's
H  = randi(H_range, n_gen, 1);   % inertia constants

D_range = [1 6];  % range of damping
Dg = randi(D_range, n_gen, 1);   % damping

Mass_gen = 2*H/omega_r;
M = diag(Mass_gen);     % mass matrix for the generators
Dg = diag(Dg);          % damping matrix for the generators



%% System Matrices
%%%% state update matrix                        - 2*n_gen x 2*n_gen
A_bar = [zeros(n_gen)                   eye(n_gen);
        -M\(L_gg-(L_gl*(L_ll\L_lg)))    -M\Dg];
%%%% input matrix                               - 2*n_gen x n_gen+n_bus
%        Pg           Pd
B_bar = [zeros(n_gen) zeros(n_gen, n_load);
        inv(M)       -M\(L_gl/L_ll)];
%%%% output matrix                              - 
%%% in terms of x_bar = [delta; omega];
C_omega       = [zeros(n_gen)     eye(n_gen)];
C_theta       = [-L_ll\L_lg       zeros(n_load,n_gen)];
C_delta_slack = zeros(1,n_states);
C_delta_slack(slack_mask) = 1;
% theta will not be an output
C_Pnet  = Lap*[eye(n_gen) zeros(n_gen);-inv(L_ll)*L_lg zeros(n_load,n_gen)]; % [L_gg - L_gl*(L_ll\L_lg)   zeros(n_gen)];  %  % [-(L_ll\L_lg)   zeros(n_load,n_gen)];
% Does not consider the power injected from generators
C_bar   = [C_omega; C_theta; C_Pnet];
C_obsv  = [C_omega;C_delta_slack; C_Pnet];

%%%% feedforward matrix
D_omega = zeros(n_gen,n_input);
D_theta = zeros(n_input);
D_delta_slack = zeros(1,n_input);
D_Pnet  = Lap*[zeros(n_gen,n_input);zeros(n_load,n_gen) -inv(L_ll)]; % [zeros(n_gen) -L_gl/(L_ll)];% %-ones(n_load);
D_bar   = [D_omega; D_theta; D_Pnet];
D_obsv  = [D_omega;D_delta_slack; D_Pnet];

%% Descretized System Model
T_sample = 0.01;
[A_bar_d, B_bar_d] = discretize_linear_model(A_bar,B_bar,T_sample);
C_obsv_d = C_obsv;
D_obsv_d = D_obsv;

disp('eigenvalues of linearized A')
disp(eig(A_bar_d).')

%% system matrix unit testing
%controllability and observability
disp('controllability')
disp(rank(ctrb(A_bar_d,B_bar_d))) % fully controllable with PID controller
disp('observability')
disp(rank(obsv(A_bar_d,C_obsv_d))) % fully observable

%% Observer Dynamics
%%% Pole Placement
P = .1 + .4*rand(2*n_gen,1);
L_obsv = place(A_bar_d.',C_obsv_d.',P).';
disp('discrete observer (A-L*C) eigenvalues: negative?')
disp(eig(A_bar_d-L_obsv*C_obsv_d).')
 
%% Simulation Initialization
x0          = zeros(n_states,1);
x0_hat      = zeros(n_states,1);
load_buses  = ones(n_bus-n_gen,1); %[zeros(n_gen,1); ones(n_bus-n_gen,1)];

%% Resilient Observer Parameters

% parameters
n       = n_states;           % # of states
n_meas  = n_gen + n_bus + 1; %n_gen +  n_gen + 1; % # of measurements
m       = n_meas;           
l       = n_input;            % # of inputs

T       = round(2*n);     % receeding horizon

N_samples      = 800; % The total number of samples to run
T_final        = N_samples*T_sample;  % Total time for simulation
% T_start_opt    = 1.5*T*T_sample; % start time for the optimization routines. Wait to collect enoguh data before calling the optimizers
tau     = 0.1;              % auxiliary model (1-reliability level)

% tapped delay
U0      = zeros(l,T);
Y0      = zeros(m,T);


[PhiT,HT,Theta_T,G_T,r_tau] = opti_params(A_bar_d,B_bar_d,C_obsv_d,T,tau);

% solver initial condition
z_0         = zeros(n+m*T,1);

%% Attack Parameters
T_start_attack = .2*T_final;  % Time to begin attack. Neede to sshow system responses with/without attacks in the same simulation
n_attack =  round(0.3*n_input);
% max_attack = 100; % maximum allowable attack per channel

%% Bad Data Detection
BDD_thresh = 0.5e-1;  % Bad data detection tolerance
[U,~,~] = svd(C_obsv_d);
U2 = U(:,n_states+1:end);


%% Auxiliary model

% % NYISO Auxiliary Models (BUILD OUT LATER)
% data_GPRs_struct = load('AuxiliaryData\GPRs_2018JantoMar.mat');
% GPRs_fieldnames = fieldnames(data_GPRs_struct);
% gprs      = data_GPRs_struct.(GPRs_fieldnames{1});
% z_bus_all = data_GPRs_struct.(GPRs_fieldnames{2});
% 
% % Measurements Data
% data_y_struct = load('AuxiliaryData\measurements_2018JantoMar.mat');
% y_meas_all = data_y_struct.measurements_nozero;
% 
% % Using the 1st 14 smallest signals
% y_means = mean(y_meas_all,1);
% N = size(y_meas_all,1);
% [~,I_sort] = sort(abs(y_means));



% Generating surrogate auxiliary model
%  Steps:
%     1. Turn off attack and turn off resilient optimizers
%     2. Run modle without attacks and collect datat
%     3. Estimate the std for each channel
%     4. Turn attack and optimizers back on and set the right parameters

n_stds = 3;  % number of standard deviations to locate auxiliary mean

% Initial Sigmas
sigma_inv_k = 1e4*(3 + 2*rand(m,1)); % surrogate inverse covariance values. Assume diagonal covariance matrix.
Sigma_inv_k = diag(sigma_inv_k);

U_y_1 = 0;

% turn off attacks and optimizers
T_start_opt = T_final*2;  % optimizers will never start
max_attack  = 0;  % zero attack

% run simulation without attacks
n_comp       = 3;% number of compopnents to retain
U_y_1 = zeros(n_meas,n_comp);
sim_out = sim('System_model_discrete');
y_obsv       = sim_out.logsout.getElement('y_obsv').Values.Data.'; % measured y for observer (y_obsv = y-D*u)
[U_y,S_y,V_y]= svd(y_obsv);
U_y_1(:,:)   = U_y(:,1:n_comp);  % dimesionality reduction
weights      = U_y_1.'*(y_obsv);
y_obsv_hat   = U_y_1*weights;
error        = y_obsv - y_obsv_hat;

error_norms = max(sum(error.^2));

% update Sigmas
sigma_inv_k(:)   = (1/error_norms)*ones(n_meas,1);%(error_norms.^(-1));
Sigma_inv_k(:,:) = diag(sigma_inv_k);

r_tau = max(sum(Sigma_inv_k*(error.^2)));

% calculate stds for y_obsv
y_stds = std(y_obsv,0,2);

% update Sigmas
sigma_inv_k(:)   = (y_stds.^(-2));
Sigma_inv_k(:,:) = diag(sigma_inv_k);


% % turn on attacks and optimizers
T_start_opt(:)    = 1.5*T*T_sample; % start time for the optimization routines. Wait to collect enoguh data before calling the optimizers
max_attack(:)     = 100;



