%% Moving Horizon Estimation TClab
% Johannes Stockhammer
% The CasADi/IPOPT set-up (nlpsol options, solver call) follows the pattern of
% the public CasADi MPC/MHE examples by M. W. Mehrez; see THIRD_PARTY.md.
clear all;
clc;
%% Import data & casadi
addpath('casadi-windows-matlabR2016a-v3.5.5');  % add casadi file path
import casadi.*                                 % import casadi
data = importdata('MHE_data.mat');              % import simulation data
est_horizon = 250;                             % estimation horizon *seconds*
%% Model parameters
t_period = 1;                                   % Discrete time period *1 seconds*
n_states = 2;                                   % number of states
n_controls = 3;                                 % number of control inputs
n_measurements = 2;                             % number of measurements
Q = diag([0.5,0.5]);                            % variance process
Qs = sqrt(Q);                                   % standard deviation process
Qsi = inv(Qs);                                  % inverse of std. dev.
R = diag([4,4]);                                % variance process
Rs = sqrt(R);                                   % standard deviation process
Rsi = inv(Rs);                                  % inverse of std. dev.
A = 1e-3;                                       % surface area not between heaters m^2
As = 2e-4;                                      % surface area between heaters m^2
m   = 0.004;                                    % mass (m1=m2) kg
UU   = 2.3;                                     % overall heat transfer coefficient W/(m^2*K)
E   = 0.9;                                      % emissivity unit 1
Sig = 5.67e-8;                                  % stefan boltzmann constant W/(m^2*K^4)
%% Initialize Variables
cp_est = zeros(1,0);                                                    
alpha1_est = zeros(1,0);                                                        
alpha2_est = zeros(1,0);
x_est = zeros(n_states,0);
status = 0;
%% Model constants
alpha1  = 0.01;                                % heat factor W/%heater
alpha1_ = SX.sym('alpha1_');
alpha2_ = SX.sym('alpha2_');
alpha2  = 0.0075;                              % heat factor W/%heater
%estimated* cp  = 500;                          % heat capacity (cp1=cp2) J/(Kg*K)
cp_ = SX.sym('cp_');                            % estimation symbol for cp
%% Define symbols for nonlinear programming problem
MEASUREMENTS = SX.sym('MEASUREMENTS',n_measurements,est_horizon);   % symbol for measurements
CONTROLS = SX.sym('CONTROLS',n_controls,est_horizon);               % symbol for control inputs
STATES = SX.sym('STATES',n_states,est_horizon);                     % symbol for states
CP_EST = SX.sym('CP_EST');                                          % symbol for CP
%% Define casadi model function
T1 = SX.sym('T1');                              % symbol for T1
T2 = SX.sym('T2');                              % symbol for T2
states_f = [T1;T2];                             % function state vecter
Q1 = SX.sym('Q1');                              % symbol for Q1
Q2 = SX.sym('Q2');                              % symbol for Q2
Tu = SX.sym('Tu');                              % symbol for Tu
controls_f = [Q1;Q2;Tu];                        % function control vector
c1  = UU*A/(m*cp_);                             % model constant 1
c2  = E*Sig*A/(m*cp_);                          % model constant 2
c3 = 1/(m*cp_);                                 % model constant 3
dT = [c1*(Tu-T1)+c2*(Tu^4-T1^4)+c3*UU*As*(T2-T1)+c3*E*Sig*A*(T2^4-T1^4)+c3*alpha1_*Q1; ...
    c1*(Tu-T2)+c2*(Tu^4-T2^4)-c3*UU*As*(T2-T1)-c3*E*Sig*A*(T2^4-T1^4)+c3*alpha2_*Q2];        % model equasion
ode = Function('ode',{states_f,controls_f,cp_,alpha1_,alpha2_},{dT});                          % casadi model function
%% Define cost function
cost_function = 0;                                                      % no initial costs
mes_init = MEASUREMENTS(:,1);                                           % initial heater temperature
state_init = STATES(:,1);                                               % initial estimated state value
g = [ mes_init-state_init ];                                            % initial constraint = init heater temperature
for j = 2:est_horizon
    mes = MEASUREMENTS(:,j);                                            % measurement at time-point j
    state = STATES(:,j);                                                % state at time-point j
    cost_function = cost_function + (mes-state)'*(mes-state)*10+0.001*(CP_EST-500)^2;           % penulty measurement-state deviation
end                                                                     % end for
%% Define constraint vector
for j = 1:est_horizon-1
    control_k = CONTROLS(:,j);                                          % control input at time-point j
    state_k = STATES(:,j);                                              % state at time-point j
    state_kp1 = STATES(:,j+1);                                          % state at time-point j+1
    state_kp1_pred = state_k+t_period*ode(state_k,control_k,CP_EST,alpha1,alpha2);    % predicted state based on model function
    g = [g; state_kp1-state_kp1_pred];                                  % hard contraint on state trajectory-model deivation
end                                                                     % end for
%% Sum up nonlinear programming problem
PARAMETERS = [reshape(MEASUREMENTS,1,numel(MEASUREMENTS)),reshape(CONTROLS,1,numel(CONTROLS))]; % input parameters for NLPP (matrix would be also allowed)
estimation_variables = [CP_EST,reshape(STATES,1,numel(STATES))];          % optimization variables (only vector allowed)
nlpp = struct('f', cost_function, 'x', estimation_variables, 'g', g, 'p', PARAMETERS);          % store nlpp in struct
%% Define solver & options
options = struct;                                 % options as struct
options.ipopt.max_iter = 100;                     % max. numerical iterations
options.ipopt.print_level = 0;
options.print_time = 0;
options.ipopt.acceptable_tol = 1e-9;
options.ipopt.acceptable_obj_change_tol = 1e-9;
solver = nlpsol('solver','ipopt',nlpp,options);   % use interior point optimizer
%% Set constraint
lbx(1) = 100;                   % lower bound on CP_EST
ubx(1) = 2000;                  % upper bound on CP_EST
lbx(1:1+numel(STATES)) = 290;   % lower bound on states
ubx(1:1+numel(STATES)) = 2000;  % upper bound on states
lbg = 0;                        % lower bound g
ubg = 0;                        % upper bound g
%% Moving Horizon Estimation
while length(data) > 4450 %est_horizon
    %% Preparing data input - move horizon
    T1m = data(2,1:est_horizon);                                % extract measurement T1
    T2m = data(3,1:est_horizon);                                % extract measurement T2
    Q1 = data(4,1:est_horizon);                                 % extract control input Q1
    Q2 = data(5,1:est_horizon);                                 % extract control input Q2
    Tu = data(6,1:est_horizon);                                 % extract control input Tu
    Y_ = [T1m;T2m];                                             % define measurement vector
    U_ = [Q1;Q2;Tu];                                            % define control input vector
    %% Run solver & extract solution
    p = [reshape(Y_,1,numel(Y_)),reshape(U_,1,numel(U_))];                                  % set up solver input parameters
    x0(1) = 500;                                                                            % initialize CP_EST
    x0(2:1+numel(STATES)) = 330;                                                            % initialize STATES
    solution = solver('x0', x0, 'lbx', lbx, 'ubx', ubx, 'lbg', lbg, 'ubg', ubg, 'p', p);    % run solver
    cp_est(1,end+1) = full(solution.x(1));                                                  % extract CP_EST solution
    est = reshape(full(solution.x(4:end)),2,[]);
    x_est(:,end+1) = est(:,end);                                                            % extract STATES solution
    %% Moving Horizon
    data(:,1) = [];                                                                         % delete 1st data entry   
    cp_est(1,end)
end
%% Plots
timespan = est_horizon:1:est_horizon+length(cp_est)-1;
value_mean = zeros(1,length(cp_est(1,:)));
%% Graph T1
subplot(3,1,1);
data = importdata('MHE_data.mat');                                          % import simulation data
plot_data = [data(2,est_horizon:1:est_horizon+length(cp_est)-1);data(3,est_horizon:1:est_horizon+length(cp_est)-1); x_est];% x_est(1,:)];
plot(timespan,plot_data);
% %% Graph T2
% subplot(4,1,2);
% data = importdata('MHE_data.mat');                                          % import simulation data
% plot_data = data(3,est_horizon:1:est_horizon+length(cp_est)-1);% x_est(2,:)];
% plot(timespan,plot_data);
%% Graph cp_eat and mean
subplot(3,1,2);
for i = 1:length(cp_est(1,:))
    value_mean(1,i) = mean(cp_est(1,:));
end
plot_data = [cp_est(1,:); value_mean(1,:)];
plot(timespan,plot_data);
%% Graph control
subplot(3,1,3);
data = importdata('MHE_data.mat');                                          % import simulation data
plot_data = [data(4,est_horizon:1:est_horizon+length(cp_est)-1); data(5,est_horizon:1:est_horizon+length(cp_est)-1)];
plot(timespan,plot_data);