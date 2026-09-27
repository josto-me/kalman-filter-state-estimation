% SPDX-FileCopyrightText: Johannes Stockhammer
% SPDX-License-Identifier: Apache-2.0
%
% ukf_2heater_T1.m  -  Unscented Kalman filter, dual heater, estimates T1 and T2 from the T1 measurement.
% The filter keeps its state in persistent variables: call it once per sample time
% (1 s); reset it with "clear ukf_2heater_T1".
%
function x_est = ukf_2heater_T1(y_k,u_k)
% Link to the source of UKF
% Course State Estimation WS2020
%% Function inputs
% u_km1 = [Q1;Q2;Tu]
% Q1 control input heater 1 %
% Q2 control input heater 2 %
% Tu ambient temperature K
% y_k = T1
% T1 sensor temperature heater 1 K
% T2 not measured
y_k = y_k(1);
u_k = reshape(u_k,[],1);
%% Function outputs
% x_est = [T1;T2]
% T1 temperature heater 1 K
% T2 temperature heater 2 K

ndimx = 2;                                  % number of states
persistent tperiod Ka Q G R x_hat_plus_km1 P_plus_km1 alpha1 alpha2 A As U E Sig c1 c2 c3
if(isempty(tperiod))
    %% Parameters
    tperiod = 1;                            % time period
    alpha1  = 0.01;                         % heat factor W/%heater
    alpha2  = 0.0075;                       % heat factor W/%heater
    cp = 500;                               % heat capacity (cp1=cp2) J/(Kg*K)
    A = 0.001;                              % surface area not between heaters m^2
    As = 0.0002;                            % surface area between heaters m^2
    m   = 0.004;                            % mass (m1=m2) kg
    U   = 10;                               % overall heat transfer coefficient W/(m^2*K)
    E   = 0.9;                              % emissivity (unit 1)
    Sig = 5.67e-8;                          % stefan boltzmann constant W/(m^2*K^4)
    c1  = U*A/(m*cp);                       % model constant 1
    c2  = E*Sig*A/(m*cp);                   % model constant 2
    c3 = 1/(m*cp);                          % model constant 3
    %% Parameter for sigma-point prediction
    Ka = 3-ndimx;                           % kappa
    %% Covariance marices
    Q = diag([0.5,0.5]);                    % variance process
    G = diag([1,1]);                        % variance process
    R = 4;                                  % variance measurement
    %% Initial values state & covariance
    x_hat_plus_km1 = [y_k;y_k];             % initial state
    P_plus_km1 = Q;                         % initial variance
    x_est = reshape(x_hat_plus_km1,1,[]);   % initial function output
else
    %% Model equations
    x_dot = @(x,u) [c1*(u(3)-x(1))+c2*(u(3)^4-x(1)^4)+c3*U*As*(x(2)-x(1))+c3*E*Sig*A*(x(2)^4-x(1)^4)+c3*alpha1*u(1); ...
                    c1*(u(3)-x(2))+c2*(u(3)^4-x(2)^4)-(c3*U*As*(x(2)-x(1)))-(c3*E*Sig*A*(x(2)^4-x(1)^4))+c3*alpha2*u(2)];       % model ODE
    h_k = @(x) x(1,:);                                                                                                      % measurement equation
    %% Calculate sigma-points
    x_plus_km1 = zeros(ndimx,2*ndimx+1);                                                                 % initialize dimension
    cholP_plus_km1 = chol(P_plus_km1);                                                                   % Cholesky decomposition (square root)
    x_plus_km1(:,1) = x_hat_plus_km1;                                                                    % 1. Sigma-Punkt
    for i = 1:ndimx
        x_plus_km1(:,i+1) = x_hat_plus_km1+sqrt(ndimx+Ka)*cholP_plus_km1(i,:)';                          % calculate sigma-points 1 to ndimx
    end                                                                                                  % end for
    for i = (ndimx+1):(2*ndimx)
        x_plus_km1(:,i+1) = x_hat_plus_km1-sqrt(ndimx+Ka)*cholP_plus_km1(i-ndimx,:)';                    % calculate sigma-points ndimx+1 to 2*ndimx
    end                                                                                                  % end for
    %% Predict sigma-points with euler
    x_minus_k = zeros(ndimx,2*ndimx+1);                                                                  % initialize dimension
    for i = 0:2*ndimx
        x_minus_k(:,i+1) = x_plus_km1(:,i+1)+x_dot(x_plus_km1(:,i+1),u_k)*tperiod;                       % predict sigma-points ndimx+1 to 2*ndimx
    end                                                                                                  % end for
    %% Calculate sigma-point weights
    w = zeros(1,2*ndimx+1);                                                                              % initialize dimension
    w(1) = Ka/(ndimx+Ka);                                                                                % calculate w0
    for i = 1:2*ndimx
        w(i+1) = 1/(2*(ndimx+Ka));                                                                       % calculate weights wi 1 to 2*ndimx
    end                                                                                                  % end for
    %% Calculate expectancy values
    x_hat_minus_k = (w*x_minus_k')';                                                                     % predict state expectancy value
    y_hat_minus_k = (w*h_k(x_minus_k)')';                                                                % predict measurement expectancy value
    %% Calculate covariance
    Pxy = zeros(2,1);                                                                                    % initialize covariance Pxy
    for i = 0:2*ndimx
        Pxy = Pxy+(x_minus_k(:,i+1)-x_hat_minus_k)*(h_k(x_minus_k(:,i+1))-y_hat_minus_k)'*w(i+1);        % calculate Pxy
    end                                                                                                  % end for
    Pyy = R;                                                                                             % initialize covariance Pyy
    for i = 0:2*ndimx
        Pyy = Pyy+(h_k(x_minus_k(:,i+1))-y_hat_minus_k)*(h_k(x_minus_k(:,i+1))-y_hat_minus_k)'*w(i+1);   % calculate Pyy
    end                                                                                                  % end for
    P_minus_k = G*Q*G';                                                                                  % initialize covariance P_minus_k
    for i = 0:2*ndimx
        P_minus_k = P_minus_k+(x_minus_k(:,i+1)-x_hat_minus_k)*(x_minus_k(:,i+1)-x_hat_minus_k)'*w(i+1); % calculate P_minus_k
    end                                                                                                  % end for
    %% Update
    K = Pxy/Pyy;                                                                                         % calculate Kalman-Gain
    x_hat_plus_k = x_hat_minus_k+K*(y_k-y_hat_minus_k);                                                  % state correction
    P_plus_k = P_minus_k-K*Pyy*K';                                                                       % covariance correction
    %% Variables for next step
    x_hat_plus_km1 = x_hat_plus_k;                                                                       % state
    P_plus_km1 = P_plus_k;                                                                               % covariance
    %% Output
    x_est = reshape(x_hat_plus_k,1,[]);                                                                  % reshape output vector
end                                                                                                      % end if
end                                                                                                      % end function
