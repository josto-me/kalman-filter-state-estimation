# Third-party software, hardware and sources

Nothing in this list is part of the repository except the adapted APMonitor code noted
below. Licenses as stated in the license files of the repositories.

## Needed to run

| Dependency | Used for | License | Source |
|---|---|---|---|
| MATLAB and Simulink (R2021a or newer) | all scripts and models | proprietary (The MathWorks, Inc.) | https://www.mathworks.com |
| MATLAB Support Package for Arduino Hardware | `Heater_Lab_IO.m` (`arduino`, `readVoltage`, `writePWMDutyCycle`) | proprietary, free add-on (The MathWorks, Inc.) | MATLAB Add-On Explorer |
| CasADi 3.5.5 for MATLAB | MHE (`code/mhe_2heater.m`), NLP solved with IPOPT | LGPL-3.0; the binary release bundles solvers such as IPOPT under their own licenses (IPOPT: Eclipse Public License) | https://github.com/casadi/casadi, https://web.casadi.org |
| Temperature Control Lab (TCLab) hardware | tasks with the real system | hardware kit by APMonitor | https://apmonitor.com/pdc/index.php/Main/ArduinoTemperatureControl |

## Adapted code

| Source | License | Where | What |
|---|---|---|---|
| APMonitor/arduino, `0_Test_Device/MATLAB/tclab.m` | Apache-2.0 | `code/Heater_Lab_IO.m`, `setupImpl` | Arduino connection with COM-port prompt, TMP36 voltage-to-temperature conversion, heater PWM functions. Wrapped in a `matlab.System` class, output in K. https://github.com/APMonitor/arduino |

## Models, methods and patterns taken from others (no code copied)

- **TCLab heater model and parameters** (energy balance, Table 1 of the report): APMonitor,
  https://apmonitor.com/pdc/index.php/Main/ArduinoModeling and
  https://apmonitor.com/pdc/index.php/Main/ArduinoModeling2 . Cited in the report and in the
  model function headers.
- **Kalman filter, EKF and UKF equations**: lecture notes of the course *State Estimation*
  (WS 2020/21), Laboratory for Systems Theory and Automatic Control, OVGU Magdeburg. Cited in
  the report ([3]-[5]) and in the filter function headers. The filter code is the author's
  own implementation.
- **CasADi/IPOPT set-up in `mhe_2heater.m`**: the options block and solver call follow the
  pattern of the public CasADi MPC/MHE examples by M. W. Mehrez
  (https://github.com/MMehrez/MPC-and-MHE-implementation-in-MATLAB-using-Casadi, no license
  file in that repository). Only the generic CasADi API usage is similar; model, cost
  function, horizon handling and data handling are the author's own.
- **Result plots in the report** were exported with matlab2tikz (BSD-2-Clause,
  https://github.com/matlab2tikz/matlab2tikz).

## Related, not used

- TCLab Python package, https://github.com/jckantor/TCLab (Apache-2.0).
