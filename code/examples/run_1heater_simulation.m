% SPDX-FileCopyrightText: Johannes Stockhammer
% SPDX-License-Identifier: Apache-2.0
%
% run_1heater_simulation.m  -  task 1 of the report without Simulink:
% the single-heater model is simulated with sensor noise, and SKF, EKF and UKF
% estimate the heater temperature T1 from the noisy measurement.

addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'filters'));
clear skf_1heater ekf_1heater ukf_1heater ode_1heater    % reset persistent filter states

%% Settings
dt  = 1;                                    % sample time s (same as the filters)
N   = 1800;                                 % number of samples (30 min)
Tu  = 296.15;                               % ambient temperature K
sigma_meas = 2;                             % sensor noise standard deviation K (R = 4)
t   = (0:N-1)*dt;                           % time s
Q1  = 100*(mod(t,600) < 300);               % heater power %: 5 min on, 5 min off
rng(1);                                     % reproducible noise

%% Simulation
x = Tu;                                     % model state T1 K, starts at ambient temperature
T_model = zeros(1,N); T_meas = zeros(1,N);
T_skf = zeros(1,N); T_ekf = zeros(1,N); T_ukf = zeros(1,N);
for k = 1:N
    u = [Q1(k); Tu];                        % input vector [Q1; Tu]
    T_model(k) = x;
    T_meas(k)  = x + sigma_meas*randn;      % noisy sensor value
    T_skf(k)   = skf_1heater(T_meas(k), u);
    T_ekf(k)   = ekf_1heater(T_meas(k), u);
    T_ukf(k)   = ukf_1heater(T_meas(k), u);
    x = x + ode_1heater(x, u)*dt;           % Euler step of the model
end

%% Plot
C = 273.15;                                 % K -> degC
figure('Name', 'Task 1: single heater');
subplot(3,1,1);
plot(t, T_meas-C, '.', t, T_model-C, 'k'); grid on;
ylabel('T_1 in °C'); legend('measurement', 'model', 'Location', 'best');
subplot(3,1,2);
plot(t, T_skf-C, t, T_ekf-C, t, T_ukf-C, t, T_model-C, 'k--'); grid on;
ylabel('T_1 in °C'); legend('SKF', 'EKF', 'UKF', 'model', 'Location', 'best');
subplot(3,1,3);
plot(t, Q1); grid on; ylabel('Q_1 in %'); xlabel('t in s');

%% Error
fprintf('RMS error to the model state: measurement %.2f K, SKF %.2f K, EKF %.2f K, UKF %.2f K\n', ...
    rms(T_meas-T_model), rms(T_skf-T_model), rms(T_ekf-T_model), rms(T_ukf-T_model));
