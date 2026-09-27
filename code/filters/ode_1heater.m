% SPDX-FileCopyrightText: Johannes Stockhammer
% SPDX-License-Identifier: Apache-2.0
%
% ode_1heater.m  -  Single-heater model dT1/dt = f(T1, Q1, Tu) (APMonitor model).
% The filter keeps its state in persistent variables: call it once per sample time
% (1 s); reset it with "clear ode_1heater".
%
function x_dot = ode_1heater(x,u)
%% 1 heater model
% Link to the source of the model
% https://apmonitor.com/pdc/index.php/Main/ArduinoModeling
%% Function inputs
% x = [T1]
% T1 heater 1 temperature K
% u = [Q1,Tu]
% Q1 heater 1 control input %
% Tu ambient temperature K
%% Function outputs
% x_dot = dT1dt
% dT1dt derivative T1 to t
%% Functions outputs
x_dot = 0;
% x_dot = [dtdT1]
% dtdT1 heater 1 derivative K/s
%% Model paramters
persistent c1 c2 c3
if(isempty(c1))
    U   = 10;                       % W/(m^2*K)
    m   = 0.004;                    % kg
    cp  = 500;                      % J/(Kg*K)
    E   = 0.9;                      % Emissivity (unit 1)
    A   = 1.2e-3;                   % m^2
    Sig = 5.67e-8;                  % W/(m^2*K^4)
    al  = 0.01;                     % W/%
    c1  = U*A/(m*cp);               % model constant 1
    c2  = E*Sig*A/(m*cp);           % model constant 2 
    c3  = al/(m*cp);                % model constant 3
end                                 % end for
dtdT1 = c1*u(2)-c1*x(1)+c2*u(2)^4-c2*x(1)^4+c3*u(1);
%% Output
x_dot = dtdT1;                      % output
end                                 % end function
