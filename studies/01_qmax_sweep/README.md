# 01 — Qmax sweep and concept demonstration

This first paper stage demonstrates the proposed energy-management concept:
changing the allowable battery-discharge budget `Qmax` produces different
fuel-cell and battery trajectories for the same step demand.

The equilibrium demand increases by 10 kW at 2.0 s and by another 10 kW at
2.2 s. The study uses a 50 ms control period, uniform N=10 prediction over
0.5 s, 100 DDP iterations, and `Qmax={72,36,18,3.6} As`.

```matlab
addpath('simulation/experiments/studies/01_qmax_sweep');
run_study
```

Raw runs are written to `results/01_qmax_sweep/ddp/`; `summary.csv`, the run
index, and the paper-facing trajectory, discharge-budget, hydrogen, and tracking
figures are written to `results/01_qmax_sweep/`.

The shared step-demand definition is `core/step_demand_scenario.m`. This stage
has no runtime dependency on `simulation/ronbun_iwai/`.
