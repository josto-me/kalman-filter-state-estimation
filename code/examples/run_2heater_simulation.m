% SPDX-FileCopyrightText: Johannes Stockhammer
% SPDX-License-Identifier: Apache-2.0
%
% run_2heater_simulation.m  -  task 2 of the report without Simulink:
% the dual-heater model is simulated with sensor noise. EKF and UKF estimate T1 and
% T2 from the T1 measurement only; the second UKF uses both measurements.

addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'filters'));
clear ekf_2heater_T1 ukf_2heater_T1 ukf_2heater_T1T2 ode_2heater    % reset persistent states

%% Settings
dt  = 1;                                    % sample time s
N   = 2400;                                 % number of samples (40 min)
Tu  = 296.15;                               % ambient temperature K
sigma_meas = 2;                             % sensor noise standard deviation K
t   = (0:N-1)*dt;                           % time s
Q1  = 100*(mod(t,800) < 400);               % heater 1 power %
Q2  = 100*(mod(t+200,800) < 400);           % heater 2 power %, shifted by 200 s
rng(1);

%% Simulation
x = [Tu; Tu];                               % model states [T1; T2] K
X_model = zeros(2,N); Y_meas = zeros(2,N);
X_ekf = zeros(2,N); X_ukf1 = zeros(2,N); X_ukf2 = zeros(2,N);
for k = 1:N
    u = [Q1(k); Q2(k); Tu];                 % input vector [Q1; Q2; Tu]
    X_model(:,k) = x;
    Y_meas(:,k)  = x + sigma_meas*randn(2,1);
    X_ekf(:,k)   = ekf_2heater_T1(Y_meas(1,k), u)';
    X_ukf1(:,k)  = ukf_2heater_T1(Y_meas(1,k), u)';
    X_ukf2(:,k)  = ukf_2heater_T1T2(Y_meas(:,k), u)';
    x = x + ode_2heater(x, u)*dt;           % Euler step of the model
end

%% Plot
C = 273.15;
figure('Name', 'Task 2: dual heater');
for i = 1:2
    subplot(3,1,i);
    plot(t, Y_meas(i,:)-C, '.', t, X_ekf(i,:)-C, t, X_ukf1(i,:)-C, t, X_ukf2(i,:)-C, t, X_model(i,:)-C, 'k--');
    grid on; ylabel(sprintf('T_%d in °C', i));
    legend('measurement', 'EKF (T1 only)', 'UKF (T1 only)', 'UKF (T1 and T2)', 'model', 'Location', 'best');
end
subplot(3,1,3);
plot(t, Q1, t, Q2); grid on; ylabel('Q in %'); xlabel('t in s'); legend('Q_1', 'Q_2');
