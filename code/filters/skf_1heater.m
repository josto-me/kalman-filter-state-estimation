% SPDX-FileCopyrightText: Johannes Stockhammer
% SPDX-License-Identifier: Apache-2.0
%
% skf_1heater.m  -  Standard Kalman filter, single heater, linearised and discretised model.
% The filter keeps its state in persistent variables: call it once per sample time
% (1 s); reset it with "clear skf_1heater".
%
function [x_k,P_plus_k]= skf_1heater(y_k,u_k)
% source of SKF: OVGU course State Estimation WS2020
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
%% Model
persistent x_hat_plus_km1 P_plus_km1 Tu0 T0 Q0 c5 c6 Ad Bd Q Cd G R u_km1
if (isempty(x_hat_plus_km1)) 
    % 1 heater model, discretised for the sample time tperiod (1 s)
    % linearized in operation point THeater0 = TAmbient0 = 296.15 K, AlphaHeater0 = 0
    U   = 10;                                                   % W/(m^2*K)
    m   = 0.004;                                                % kg
    cp  = 500;                                                  % J/(Kg*K)
    E   = 0.9;                                                  % Emissivity (Unit 1)
    A   = 1.2e-3;                                               % m^2
    Sig = 5.67e-8;                                              % W/(m^2*K^4) Stefan-Boltzmann constant
    al  = 0.01;                                                 % W/%
    tperiod = 1;                                              % timer period
    Tu0 = 296.15;                                               % linearization point room temperature K
    T0 = Tu0;                                                   % linearization point ambient temperature K
    Q0 = 0;                                                     % linearization point heater control %
    c5 = (U*A+E*Sig*A*4*Tu0^3)/(m*cp);                          % model constant 5
    c6 = al/(m*cp);                                             % model constant 6
    Ad = [ exp(-c5*tperiod) ];                                  % state matrix time-discrete
    Bd = [ (c6-c6*exp(-c5*tperiod))/c5, 1-exp(-c5*tperiod) ];   % input matrix time-discret
    Cd = [ 1 ];                                                 % output matrix
    Q = [ 0.5 ];                                                % process noise
    G = [ 1 ];                                                  
    R = [ 4 ];                                                  % measurement noise
    %% Initialisation 
    x_k = y_k;                                                  % initial function state
    P_plus_k = Q;                                               % initial function covariance
    x_hat_plus_km1 = y_k-T0;                                    % initial state, relation to working point
    P_plus_km1 = Q;                                             % initial variance = variance measurement
    u_km1 = [u_k(1)-Q0; u_k(2)-Tu0];                            % relate control input to working point
else                                                            % end if
    %% Input scaling
    y_k = y_k-T0;                                                   % relate measurement to working point
    %% Predictor
    x_hat_minus_k = Ad*x_hat_plus_km1+Bd*u_km1;                     % state equation
    y_hat_minus_k = Cd*x_hat_minus_k;                               % measurement estimation
    P_minus_k = Ad*P_plus_km1*Ad'+G*Q*G';                           % (co)variance
    %% Corrector
    K = P_minus_k*Cd'/(Cd*P_minus_k*Cd'+R);                         % calculate Kalman Gain
    x_hat_plus_k = x_hat_minus_k+K*(y_k-y_hat_minus_k);             % state correction
    P_plus_k = P_minus_k-K*Cd*P_minus_k;                            % (co)variance correction
    %% Function output
    x_k = x_hat_plus_k+T0;                                          % relate temperature to working point
    %% Variables for the next step
    x_hat_plus_km1 = x_hat_plus_k;                                  % state
    P_plus_km1 = P_plus_k;                                          % (co)variance
    u_km1 = [u_k(1)-Q0; u_k(2)-Tu0];                                  % relate control input to working point
end                                                                 % end if
end                                                                 % end function
