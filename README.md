# Kalman Filter State Estimation

[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.22984394-blue.svg)](https://doi.org/10.5281/zenodo.22984394) [![Code: Apache-2.0](https://img.shields.io/badge/code-Apache--2.0-blue.svg)](LICENSE) [![Report & data: CC BY 4.0](https://img.shields.io/badge/report%20%26%20data-CC%20BY%204.0-lightgrey.svg)](LICENSE-CC-BY-4.0.txt) [![Cite](https://img.shields.io/badge/cite-CITATION.cff-green.svg)](CITATION.cff)

Zustandsschätzung mit Kalman-Filtern und Moving Horizon Estimation am Temperature Control Lab (Studienprojekt, OVGU Magdeburg).

Course project report in *State Estimation* at the Laboratory for Systems Theory and
Automatic Control, Otto-von-Guericke University Magdeburg. Standard, Extended and Unscented
Kalman filters (SKF, EKF, UKF) and a Moving Horizon Estimator (MHE, with CasADi/IPOPT) are
implemented in MATLAB/Simulink and applied to the
[Temperature Control Lab](https://apmonitor.com/pdc/index.php/Main/ArduinoTemperatureControl)
(TCLab, APMonitor): an Arduino shield with two transistor heaters and two TMP36 temperature
sensors.

The report is in [`report/state-estimation-report.pdf`](report/state-estimation-report.pdf).

## Tasks

The report is split into three tasks:

| Task | Plant | Estimators |
|---|---|---|
| 1 | Simulated 1-heater model + noisy sensor | SKF (linearised at 23 °C, discretised), EKF, UKF, estimating T1 (`code/examples/run_1heater_simulation.m`) |
| 2 | Simulated 2-heater model + noisy sensors | EKF and UKF1 estimate T1, T2 from T1 only; UKF2 uses T1 and T2 (`code/examples/run_2heater_simulation.m`) |
| 3 | Real TCLab over USB | UKF1 / UKF2 as in task 2, run in real time next to the model |
| 3 | Real TCLab over USB | Record Q1, Q2, T1, T2, then MHE offline: estimate heat capacity cp and T1, T2 (`code/mhe_2heater.m`) |

The heater model is the energy balance published by APMonitor (convection, radiation,
heat exchange between the two heaters, heater power `alpha * Q`). In the report, the filters are MATLAB
Function blocks inside Simulink models; `code/filters/` holds their code as plain MATLAB
functions. The hardware connection is a MATLAB System block
([`code/Heater_Lab_IO.m`](code/Heater_Lab_IO.m)) built on the MATLAB Support Package for Arduino
Hardware.

The MHE is formulated as one nonlinear program per horizon: states on a 250 s horizon
(1 s steps, explicit Euler) plus the unknown cp as decision variables, the model as equality
constraints, and a least-squares cost on the measurement residuals with a small penalty
on the deviation of cp from its nominal value of 500 J/(kg K). IPOPT solves it, then the
horizon moves on by one sample.

## Results (from the report)

- **Task 1 (simulation):** all three filters reduce the sensor noise by a similar amount.
  The SKF on the linearised model does as well as EKF and UKF here, because the 1-heater
  model is only weakly nonlinear around the working point.
- **Task 2 (simulation):** T2 estimated from the T1 measurement alone follows the simulated
  T2 closely. UKF2 (both sensors) is no better than UKF1, since plant and filters use the
  same model.
- **Task 3 (real TCLab):** there is a clear model mismatch. After the heater is switched off
  the measured temperature keeps rising for a while, which the model does not reproduce
  (no heat capacity or heat transfer of transistor, sensor and heat sink, no heat loss
  through the pins to the board, heat factor alpha not exact). The UKF2 estimate of T2
  deviates from the measured T2.
- **MHE:** the estimated heat capacity oscillates around roughly 650 J/(kg K) instead of the
  nominal 500 J/(kg K). The report notes that fitting cp alone may not be the right way to
  correct the model, and names correcting alpha1/alpha2 and the model structure as
  next steps.

## Contents

```
report/state-estimation-report.pdf   the report (PDF)
code/filters/ode_1heater.m              single-heater model
code/filters/skf_1heater.m              standard Kalman filter, single heater (linearised model)
code/filters/ekf_1heater.m              extended Kalman filter, single heater
code/filters/ukf_1heater.m              unscented Kalman filter, single heater
code/filters/ode_2heater.m              dual-heater model
code/filters/ekf_2heater_T1.m           EKF, dual heater, T1 and T2 from the T1 measurement
code/filters/ukf_2heater_T1.m           UKF, dual heater, T1 and T2 from the T1 measurement
code/filters/ukf_2heater_T1T2.m         UKF, dual heater, T1 and T2 from both measurements
code/examples/run_1heater_simulation.m  task 1 without Simulink: model + SKF/EKF/UKF, plots
code/examples/run_2heater_simulation.m  task 2 without Simulink: model + EKF/UKF, plots
code/Heater_Lab_IO.m                    MATLAB System block for the TCLab (Arduino I/O)
code/mhe_2heater.m                      MHE with CasADi, estimates cp from MHE_data.mat
data/MHE_data.mat                       recorded TCLab run used by mhe_2heater.m
```

`data/MHE_data.mat` (MAT-file v7.3/HDF5) holds one 6 x N double matrix `MHE_data`, one column
per second: row 1 time in s, rows 2-3 measured T1, T2 in K, rows 4-5 heater inputs Q1, Q2
in %, row 6 ambient temperature (constant 296.15 K).

## How to run

Requirements:

- MATLAB (the filters and examples use no toolbox; `rms` needs the Signal Processing
  Toolbox or can be replaced by `sqrt(mean(x.^2))`).
- For the hardware block: Simulink, the *MATLAB Support Package for Arduino Hardware* and a
  TCLab (Arduino with the TCLab shield and its 5 V heater supply). `Heater_Lab_IO.m`
  connects with board type `'Uno'`.
- For the MHE: [CasADi](https://web.casadi.org/) 3.5.5 for MATLAB (the script adds the
  folder `casadi-windows-matlabR2016a-v3.5.5` next to it to the path; adapt the `addpath`
  line for other versions or operating systems).

Filters on the simulated heaters: run `code/examples/run_1heater_simulation.m` or
`code/examples/run_2heater_simulation.m`. Each filter keeps its state in persistent
variables and is called once per second; the scripts reset them with `clear`.

MHE on the recorded data:

1. Put the `code/` and `data/` folders on the MATLAB path, e.g.
   `addpath('code','data')` from the repository root.
2. Run `code/mhe_2heater.m`. It solves one NLP per second of data, which takes a
   while, and plots T1/T2, the estimated cp and the heater inputs.

`Heater_Lab_IO.m` can be used as a MATLAB System block in your own Simulink model to read T1, T2
and set Q1, Q2 on a TCLab.

## Known issues and limitations

- The code in this form (plain MATLAB functions and the two example scripts) is not tested.
  Course project code, not maintained.
- The simple TCLab model does not describe the real sensor dynamics (see results), so the
  hardware estimates are limited by model mismatch rather than by the filters.
- **Known issue, `ekf_2heater_T1.m`:** the second row of the Jacobian `dfdx` is not the
  derivative of dT2/dt. Its first element should be `c3*U*As + c2*4*x(1)^3`, and the terms
  `-c1 - c2*4*x(2)^3` are missing in its second element. The report results were
  computed with this code.
- `ukf_1heater.m` also adds Q to Pxy and Pyy. The dual-heater UKFs use a Cholesky square root
  for the sigma points, as in the real-time models of the report; the simulation models of
  the report use an element-wise square root.
- The MHE runs offline on recorded data, not in the loop. The loop end is hard-coded
  (`while length(data) > 4450`) for a 5000 s recording, which gives estimates from 250 s
  to about 800 s as in Fig. 10 of the report. `MHE_data.mat` is the recording used for
  Fig. 10.

## License

- Code in `code/`: **Apache License 2.0**, see [`LICENSE`](LICENSE) and
  [`NOTICE`](NOTICE). `code/Heater_Lab_IO.m` contains code adapted from APMonitor (Apache-2.0),
  see [`THIRD_PARTY.md`](THIRD_PARTY.md).
- The report (`report/`), the data in `data/` and this README:
  **CC BY 4.0**, see [`LICENSE-CC-BY-4.0.txt`](LICENSE-CC-BY-4.0.txt).

You may use, change and share everything, also commercially. When you pass it on or
publish something based on it, credit it as:

> Johannes Stockhammer, "State Estimation with Kalman Filter and Moving Horizon Estimation on the Temperature Control Lab", course project report, Otto-von-Guericke University Magdeburg, revised edition, Zenodo (2026), https://doi.org/10.5281/zenodo.22984395

GitHub shows the same citation under "Cite this repository" (from [`CITATION.cff`](CITATION.cff)).

## Dependencies

Needed to run: MATLAB/Simulink and the MATLAB Support Package for Arduino Hardware
(proprietary, MathWorks), CasADi (LGPL-3.0), TCLab hardware (APMonitor). Details and
licenses in [`THIRD_PARTY.md`](THIRD_PARTY.md).

## Trademarks

MATLAB and Simulink are registered trademarks of The MathWorks, Inc. Arduino is a trademark
of Arduino SA. TCLab / Temperature Control Lab refers to the kit by APMonitor. The names are
used only to identify the software and hardware.

## Author

Johannes Stockhammer
