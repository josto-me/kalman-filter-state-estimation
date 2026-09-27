% SPDX-FileCopyrightText: Johannes Stockhammer
% SPDX-License-Identifier: Apache-2.0
%
% ukf_1heater.m  -  Unscented Kalman filter, single heater, nonlinear model (Euler, 1 s).
% The filter keeps its state in persistent variables: call it once per sample time
% (1 s); reset it with "clear ukf_1heater".
%
function [x_hat_plus_k,P_plus_k] = ukf_1heater(y_k,u_k)
% source of UKF: OVGU course State Estimation WS2020
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

ndimx = 1;                                              % number of states
persistent x_hat_plus_km1 P_plus_km1 c1 c2 c3 Q G R Ka tperiod u_km1
if(isempty(x_hat_plus_km1))         
    %% Parameters
    tperiod = 1;                                      % s
    U   = 10;                                           % W/(m^2*K)
    m   = 0.004;                                        % kg
    cp  = 500;                                          % J/(Kg*K)
    E   = 0.9;                                          % emissivity (unit 1)
    A   = 1.2e-3;                                       % m^2
    Sig = 5.67e-8;                                      % W/(m^2*K^4)
    al  = 0.01;                                         % W/%
    c1  = U*A/(m*cp);                                   % model constant 1
    c2  = E*Sig*A/(m*cp);                               % model constant 2
    c3  = al/(m*cp);                                    % model constant 3
    %% Parameter for sigma-point prediction
    Ka = 3-ndimx;                                       % kappa
    %% Covariance matrices
    Q =  0.5 ;                                          % variance process
    G =  1 ;                                            % variance process
    R =  4 ;                                            % variance process
    %% Initial values state & covariance
    u_km1 = u_k;                                        % initial control input
    x_hat_plus_k = y_k;                                 % initial state output
    P_plus_k = Q;                                       % initial covariance
    x_hat_plus_km1 = y_k;                               % initial state
    P_plus_km1 = Q;                                     % initial covariance P
else
    %% Model equation
    dxdt = @(x,u) c1*u(2)-c1*x+c2*u(2)^4-c2*x^4+c3*u(1);                                                    % ODE MISO
    h_k = @(x) x;                                                                                           % output equation
    %% Calculate sigma-points
    x_plus_km1 = zeros(ndimx,2*ndimx+1);                                                                    % initialize dimension
    x_plus_km1(:,1) = x_hat_plus_km1;                                                                       % 1st sigma-point
    for i = 1:ndimx
        x_plus_km1(:,i+1) = x_hat_plus_km1+sqrt(ndimx+Ka)*sqrt(P_plus_km1(i,:)');                           % calculate sigma-points 1 to ndimx
    end                                                                                                     % end for
    for i = (ndimx+1):(2*ndimx)
        x_plus_km1(:,i+1) = x_hat_plus_km1-sqrt(ndimx+Ka)*sqrt(P_plus_km1(i-ndimx,:)');                     % calculate sigma-points ndimx+1 to 2*ndimx
    end                                                                                                     % end for
    %% Predict sigma-points with euler
    x_minus_k = zeros(ndimx,2*ndimx+1);                                                                     % initialize dimension
    for i = 0:2*ndimx
        x_minus_k(:,i+1) = x_plus_km1(:,i+1)+dxdt(x_plus_km1(:,i+1),u_km1)*tperiod;                         % predict sigma-points ndimx+1 to 2*ndimx
    end                                                                                                     % end for
    %% Calculate sigma-point weights
    w = zeros(1,2*ndimx+1);                                                                                 % initialize dimension
    w(1) = Ka/(ndimx+Ka);                                                                                   % calculate w0
    for i = 1:2*ndimx
        w(i+1) = 1/(2*(ndimx+Ka));                                                                          % calculate weights wi 1 to 2*ndimx
    end                                                                                                     % end for
    %% Calculate expectancy values
    x_hat_minus_k = (w*x_minus_k')';                                                                        % predict state expectancy value
    y_hat_minus_k = (w*h_k(x_minus_k)')';                                                                   % predict measurement expectancy value
    %% Calculate covariance
    Pxy = Q;                                                                                                % initialize covariance Pxy
    for i = 0:2*ndimx
        Pxy = Pxy+(x_minus_k(:,i+1)-x_hat_minus_k)*(h_k(x_minus_k(:,i+1))-y_hat_minus_k)'*w(i+1);           % calculate Pxy
    end                                                                                                     % end for
    Pyy = R+Q;                                                                                              % initialize covariance Pyy
    for i = 0:2*ndimx
        Pyy = Pyy+(h_k(x_minus_k(:,i+1))-y_hat_minus_k)*(h_k(x_minus_k(:,i+1))-y_hat_minus_k)'*w(i+1);      % calculate Pyy
    end                                                                                                     % end for
    P_minus_k = G*Q*G';                                                                                     % initialize covariance P_minus_k
    for i = 0:2*ndimx
        P_minus_k = P_minus_k+w(i+1)*(x_minus_k(:,i+1)-x_hat_minus_k)*(x_minus_k(:,i+1)-x_hat_minus_k)';    % calculate P_minus_k
    end                                                                                                     % end for
    %% Update
    K = Pxy/Pyy;                                                                                            % calculate Kalman-Gain
    x_hat_plus_k = x_hat_minus_k+K*(y_k-y_hat_minus_k);                                                     % state correction
    P_plus_k = P_minus_k-K*Pyy*K';                                                                          % covariance correction
    %% Variables for next step
    x_hat_plus_km1 = x_hat_plus_k;                                                                          % state
    P_plus_km1 = P_plus_k;                                                                                  % covariance
    u_km1 = u_k;                                                                                            % store control input
end                                                                                                         % end if
end                                                                                                         % end function
