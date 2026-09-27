% SPDX-FileCopyrightText: Johannes Stockhammer
% SPDX-License-Identifier: Apache-2.0
%
% ode_2heater.m  -  Dual-heater model [dT1/dt; dT2/dt] = f(T1, T2, Q1, Q2, Tu) (APMonitor model).
% The filter keeps its state in persistent variables: call it once per sample time
% (1 s); reset it with "clear ode_2heater".
%
function x_dot = ode_2heater(x,u)
% Link to the source of the model
% https://apmonitor.com/pdc/index.php/Main/ArduinoModeling2
%% Function inputs
% x = [T1;T2]
% u = [Q1;Q2;Tu]
% Q1 control input heater 1,% (1 W max.)
% Q2 control input heater 2 % (0.75 W max.)
% T1 temperature heater 1 K
% T2 temperature heater 2 K
% Tu ambient temperature K
T1 = x(1);  % 1st state = T1
T2 = x(2);  % 2nd state = T2
Q1 = u(1);  % 1st control input = Q1
Q2 = u(2);  % 2nd control input = Q2
Tu = u(3);  % 3rd control input = Tu
%% Function outputs
% dT1dt derivative of T1
% dT2dt derivative of T2
%% Model constants
persistent alpha1 alpha2 A As U E Sig c1 c2 c3                      % persistent variables
if(isempty(alpha1))
    alpha1  = 0.01;                                                 % heat factor W/%heater
    alpha2  = 0.0075;                                               % heat factor W/%heater
    cp  = 500;                                                      % heat capacity (cp1=cp2) J/(Kg*K)
    A = 0.001;                                                      % surface area not between heaters m^2
    As = 0.0002;                                                    % surface area between heaters m^2
    m   = 0.004;                                                    % mass (m1=m2) kg
    U   = 10;                                                       % overall heat transfer coefficient W/(m^2*K)
    E   = 0.9;                                                      % emissivity (unit 1)
    Sig = 5.67e-8;                                                  % stefan boltzmann constant W/(m^2*K^4)
    c1  = U*A/(m*cp);                                               % constant 1 
    c2  = E*Sig*A/(m*cp);                                           % constant 2 
    c3 = 1/(m*cp);                                                  % constant 3
end                                                                 % end if
%% Model equations
Qc12 = U*As*(T2-T1);                                                % convective heat transfer
Qr12 = E*Sig*A*(T2^4-T1^4);                                         % radiative heat transfer
x_dot = [c1*(Tu-T1)+c2*(Tu^4-T1^4)+c3*Qc12+c3*Qr12+c3*alpha1*Q1;    % 1st heater
         c1*(Tu-T2)+c2*(Tu^4-T2^4)-c3*Qc12-c3*Qr12+c3*alpha2*Q2];   % 2nd heater
end                                                                 % end function
