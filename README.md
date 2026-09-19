# FCV-MPC numerical experiments

This repository contains the MATLAB/Simulink code and generated data for the
numerical studies in `FCV_MPC.pdf`. The numbered study folders follow the
intended order of the paper.

## Paper validation flow

| Stage | Question addressed | Entry point | Results | Current status |
|---|---|---|---|---|
| `01_qmax_sweep` | Can the proposed energy-management structure shape different FC/battery trajectories by changing `Qmax`? | `studies/01_qmax_sweep/run_study.m` | `results/01_qmax_sweep/` | Four-condition concept demonstration available |
| `02_method_comparison` | How does the proposed method compare with conventional energy management under matched conditions? | `studies/02_method_comparison/run_study.m` | `results/02_method_comparison/` | Diagnostic; comparison conditions and convergence still need refinement |
| `03_udds` | How does the method behave on a driving-cycle-derived demand? | `studies/03_udds/run_study.m` | `results/03_udds/` | Representative-interval results available; final paper framing is open |

`studies/solver_validation/` is a supporting numerical check rather than a
separate stage of the paper. It compares DDP and fmincon solution quality and
writes to `results/solver_validation/`.

Run one stage by adding only its folder, for example:

```matlab
addpath('simulation/fcv-mpc-energy-management/studies/01_qmax_sweep');
run_study

rmpath('simulation/fcv-mpc-energy-management/studies/01_qmax_sweep');
addpath('simulation/fcv-mpc-energy-management/studies/02_method_comparison');
run_study

rmpath('simulation/fcv-mpc-energy-management/studies/02_method_comparison');
addpath('simulation/fcv-mpc-energy-management/studies/03_udds');
run_study
```

## Directory map

- `studies/01_qmax_sweep/`: step-demand concept demonstration.
- `studies/02_method_comparison/`: step-demand conventional-method comparison.
- `studies/03_udds/`: driving-cycle-derived validation.
- `core/`: shared scenarios, demand construction, method drivers, and metrics.
- `methods/`: DDP, direct-shooting fmincon, and low-pass implementations.
- `results/<stage>/`: raw data, summary tables, and figures matching each stage.

The control period is 50 ms and the plant RK4 step is 5 ms. Uniform N=10 and
non-uniform N=5 horizons both cover 0.5 s. Battery-use and hydrogen comparisons
are accepted only when demand, initial state, constraints, and battery use are
comparable. Constraint violations and solver failures remain visible.

MATLAB/Simulink caches such as `slprj/` and `*.slxc` are not part of the public
repository. Historical material under `simulation/archive/` and
`simulation/ronbun_iwai/` is not a runtime dependency.
