% SPDX-FileCopyrightText: Johannes Stockhammer
% SPDX-License-Identifier: Apache-2.0
%
% ekf_2heater_T1.m  -  Extended Kalman filter, dual heater, estimates T1 and T2 from the T1 measurement.
% The filter keeps its state in persistent variables: call it once per sample time
% (1 s); reset it with "clear ekf_2heater_T1".
%
function x_est = ekf_2heater_T1(y_k,u_k)
% source of EKF: OVGU course State Estimation WS2020
%% Function inputs
% y_k = [T1,T2]
% T1 heater 1 temperature K
% T2 heater 2 temperature K
% u_k = [Q1;Q2;Tu]
% Q1 heater 1 control input %
% Q2 heater 2 control input %
% Tu ambient temperature K
%y_k = reshape(y_k,[],1);
y_k = y_k(1);
u_k = reshape(u_k,[],1);
%% Function outputs
% x_hat_plus_k = [T1,T2]
% P_plus_k = [ P11,P12 ;
%              P21,P22 ];
persistent x_hat_plus_km1 P_plus_km1 tperiod Qd G R c1 c2 c3 U As alpha1 alpha2 n_states n_measurements
if(isempty(x_hat_plus_km1))
    %% Model parameters
    n_states = 2;
    n_measurements = numel(y_k);
    tperiod = 1;                            % s
    alpha1  = 0.01;                         % heat factor W/%heater
    alpha2  = 0.0075;                       % heat factor W/%heater
    A   = 0.001;                            % surface area not between heaters m^2
    As  = 0.0002;                           % surface area between heaters m^2
    m   = 0.004;                            % mass (m1=m2) kg
    U   = 10;                               % overall heat transfer coefficient W/(m^2*K)
    cp  = 500;                              % heat capacity (cp1=cp2) J/(Kg*K)
    E   = 0.9;                              % emissivity unit 1
    Sig = 5.67e-8;                          % stefan boltzmann constant W/(m^2*K^4)
    c1  = U*A/(m*cp);                       % model constant 1
    c2  = E*Sig*A/(m*cp);                   % model constant 2
    c3  = 1/(m*cp);                         % model constant 3
    %% Covariance
    Q = diag([0.5,0.5]);                    % variance process
    Qd = Q/tperiod;                         % calculate Qd
    G = diag([1,1]);       
    %R = diag([4,4])                        % variance measurement
    R = 4;
    %% Initial values
    x_hat_plus_km1 = [ y_k; y_k ];          % initial state = measured temperature
    P_plus_km1 = Qd;                        % initial state variance = measurement variance
    x_est = reshape(x_hat_plus_km1,1,[]);   % initial function output
else
    %% Model functions
    dxdt = @(x,u) [ c1*(u(3)-x(1))+c2*(u(3)^4-x(1)^4)+c3*U*As*(x(2)-x(1))+c2*(x(2)^4-x(1)^4)+c3*alpha1*u(1) ; ...
                    c1*(u(3)-x(2))+c2*(u(3)^4-x(2)^4)-c3*U*As*(x(2)-x(1))-c2*(x(2)^4-x(1)^4)+c3*alpha2*u(2) ];  % model ODE
    dfdx = @(x) [ -c1-c2*4*x(1)^3-c3*U*As-c2*4*x(1)^3    c3*U*As+c2*4*x(2)^3 ;
                  -c1-c2*4*x(1)^3+c3*U*As+c2*4*x(1)^3   -c3*U*As-c2*4*x(2)^3 ];  % Jacobi matrix
    h = @(x) [1 0]*x;                                    % output function
    dxdt_hat_km1 = dxdt(x_hat_plus_km1,u_k);                                                                    % dxdt at time km1
    %% Partial derivatives dfdx and dhdx
    A = dfdx(x_hat_plus_km1);                            % derivative of the model ODE at x_hat_km1
    dhdx = [1 0];                                        % derivative of the output equation h(x) with respect to x
    %% Covariance of the state prediction
    dPdt_km1 = A*P_plus_km1+P_plus_km1*A'+G*Qd*G';       % derivative of P_minus_k at km1
    %% Integrate with Euler
    x_hat_minus_k = x_hat_plus_km1+dxdt_hat_km1*tperiod; % integrate state
    P_minus_k = P_plus_km1+dPdt_km1*tperiod;             % integrate covariance
    %% Measurement prediction
    y_hat_k = h(x_hat_minus_k);                          % measurement prediction = predicted state
    %% Kalman Gain
    Pxy = P_minus_k*dhdx';
    Pyy = dhdx*P_minus_k*dhdx'+R; 
    K_plus_k = Pxy/Pyy;                                  % Kalman Gain
    %% State and covariance correction
    x_hat_plus_k = x_hat_minus_k+K_plus_k*(y_k-y_hat_k); % weighted state
    P_plus_k = P_minus_k-K_plus_k*Pyy*K_plus_k';         % weighted covariance
    %% Variables for the next step
    x_hat_plus_km1 = x_hat_plus_k;                       % state
    P_plus_km1 = P_plus_k;                               % covariance
    %% Output
    x_est = reshape(x_hat_plus_k,1,[]); 
end                                                      % end if
end                                                      % end function
