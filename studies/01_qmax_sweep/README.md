# 01 — Qmax sweep and concept demonstration

This first paper stage demonstrates the proposed energy-management concept:
changing the allowable battery-discharge budget `Qmax` produces different
fuel-cell and battery trajectories for the same step demand.

The equilibrium demand increases by 10 kW at 2.0 s and by another 10 kW at
2.2 s. The study uses a 50 ms control period, uniform N=10 prediction over
0.5 s, 100 DDP iterations, and `Qmax={72,36,18,3.6} As`.

```matlab
addpath('simulation/fcv-mpc-energy-management/studies/01_qmax_sweep');
run_study
```

Raw runs are written to `results/01_qmax_sweep/ddp/`; `summary.csv`, the run
index, and the paper-facing trajectory, discharge-budget, hydrogen, and tracking
figures are written to `results/01_qmax_sweep/`.

The figures are numbered in manuscript order:

1. `01_power_tracking.png`: supplied power and demand (manuscript Fig. 3)
1a. `01a_tracking_error.png`: diagnostic `Psys-Pref` companion to the first figure
2. `02_control_inputs.png`: compressor voltage, stack current, and battery current
3. `03_battery_discharge.png`: cumulative battery discharge
4. `04_hydrogen_consumption.png`: hydrogen comparison

The plotting script uses a common IEEE-paper style for canvas proportions,
Times New Roman typography, line widths, axes, colors, and legends.
`summary.csv` records full-window RMSE, maximum absolute error, mean error,
final error, and post-step RMSE. Plot styling must not hide a tracking discrepancy.

The shared step-demand definition is `core/step_demand_scenario.m`. This stage
has no runtime dependency on `simulation/ronbun_iwai/`.
