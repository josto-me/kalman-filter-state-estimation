% SPDX-FileCopyrightText: Johannes Stockhammer
% SPDX-License-Identifier: Apache-2.0
%
% ekf_1heater.m  -  Extended Kalman filter, single heater, nonlinear model (Euler, 1 s).
% The filter keeps its state in persistent variables: call it once per sample time
% (1 s); reset it with "clear ekf_1heater".
%
function [x_hat_plus_k,P_plus_k]= ekf_1heater(y_k,u_k)
% source of EKF: OVGU course State Estimation WS2020
%% Function inputs
% y_k = [T1]
% sensor T1 heater 1 temperature K
% u_km1 = [Q1,Tu]
% Q1 heater 1 control input %
% Tu ambient temperature K
%% Function outputs
% x_k = [T1]
% P_plus_k
% T1 estimated heater 1 temperature K
% P Covariance matrix estimated state
%% Constants of the differential equation
persistent x_hat_plus_km1 P_plus_km1 tperiod Qd G R c1 c2 c3 u_km1
if(isempty(x_hat_plus_km1))
    %% Model parameters
    U   = 10;                                           % W/(m^2*K)
    m   = 0.004;                                        % kg
    cp  = 500;                                          % J/(Kg*K)
    E   = 0.9;                                          % Emissivity (unit 1)
    A   = 1.2e-3;                                       % m^2
    Sig = 5.67e-8;                                      % W/(m^2*K^4)
    al  = 0.01;                                         % W/%
    tperiod = 1;                                      % s
    %% (Co)variances
    Q =  0.5;                                          % variance process
    Qd = Q/tperiod;
    G = 1;                          
    R = 4;                                              % variance process
    %% Initial values
    u_km1 = u_k;                                        % initial control input
    x_hat_plus_k = y_k;                                 % initial state output
    P_plus_k = Q;                                       % initial covariance
    x_hat_plus_km1 = y_k;                               % initial state, relation to working point
    P_plus_km1 = Q;                                     % initial variance = variance measurement
    %% Model constants
    c1  = U*A/(m*cp);                                   % model constant 1
    c2  = E*Sig*A/(m*cp);                               % model constant 2
    c3  = al/(m*cp);                                    % model constant 3
else
    %% Model equations
    dxdt = @(x,u) c1*u(2)-c1*x+c2*u(2)^4-c2*x^4+c3*u(1);    % ODE MISO
    dfdx = @(x) -c1-c2*4*x^3;                               % jacobi
    dhdx = 1;                                               % derivative von h(x) to x
    h = @(x) x;                                             % measurement equation
    %% Partial derivative dfdx
    A = dfdx(x_hat_plus_km1);                               % derivative at state x_hat_km1
    %% Covariance prediction
    dPdt_km1 = A*P_plus_km1+P_plus_km1*A'+G*Qd*G';           % derivative P_minus_k at k-1
    %% Integrate with euler
    dxdt_hat_km1 = dxdt(x_hat_plus_km1,u_km1);                % dxdt at time-point km1
    x_hat_minus_k = x_hat_plus_km1+dxdt_hat_km1*tperiod;    % use euler for state
    P_minus_k = P_plus_km1+dPdt_km1*tperiod;                % use euler for covariance
    %% Measurement prediction
    y_hat_k = h(x_hat_minus_k);                             % measurement prediction
    %% Kalman Gain
    Pxy = P_minus_k*dhdx';
    Pyy = dhdx*P_minus_k*dhdx'+R;
    K_plus_k = Pxy/Pyy;                                     % calculate Kalman-Gain
    %% State and covariance correction
    x_hat_plus_k = x_hat_minus_k+K_plus_k*(y_k-y_hat_k);    % weighted state
    P_plus_k = P_minus_k-K_plus_k*Pyy*K_plus_k';            % weighted covariance
    %% Variables for the next step
    x_hat_plus_km1 = x_hat_plus_k;                          % state
    P_plus_km1 = P_plus_k;                                  % covariance
    u_km1 = u_k;                                             % store control input
end                                                         % end if
end                                                         % end function
