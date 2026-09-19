# FCV-MPC numerical experiments

This directory contains the MATLAB/Simulink code and generated data used for
the numerical studies in `FCV_MPC.pdf`. Each paper study has the same three
entry points: a scenario file, `run_study.m`, and plotting/reporting code.

## Paper studies

| Study | Paper role | Entry point | Results |
|---|---|---|---|
| `step_demand` | Sec. V-A, Figs. 3--6: IDC behavior for four Qmax values | `studies/step_demand/run_study.m` | `results/step_demand/` |
| `udds` | Sec. V-B/C: method and computation-time comparison | `studies/udds/run_study.m` | `results/udds/` |
| `solver_validation` | DDP/fmincon solution-quality checks | `studies/solver_validation/check_solver_snapshot.m` | `results/solver_validation/` |

Run one study by adding only its folder. For example:

```matlab
addpath('simulation/experiments/studies/step_demand');
run_study

rmpath('simulation/experiments/studies/step_demand');
addpath('simulation/experiments/studies/udds');
run_study
```

Both studies call the same code in `core/` and the same implementations in
`methods/`. Scenario files change the demand, interval, Qmax, and declared
comparison settings.

## Directory map

- `studies/`: paper-facing configurations, run scripts, and figure scripts.
- `core/`: shared demand construction, method drivers, and metric evaluation.
- `methods/ddp/`: DDP solver and N=5/N=10 Simulink models.
- `methods/fmincon/`: direct-shooting reference implementation.
- `methods/lowpath/`: frequency-decoupling benchmark.
- `results/<study>/`: raw method data, tables, and figures for each study.

`slprj/` and `*.slxc` are MATLAB-generated caches and should be excluded from
the public repository. Historical development files remain outside this tree
under `simulation/archive/` and `simulation/ronbun_iwai/`; they are not runtime
dependencies of these studies.

## Shared comparison conditions

The control period is 50 ms and the plant RK4 step is 5 ms. Uniform N=10 and
non-uniform N=5 horizons both cover 0.5 s. The first non-uniform step is 50 ms.
DDP and fmincon hold the first optimized input for one control period.

Battery use is reported in As and independently checked against SOC change.
Hydrogen comparisons are accepted only when demand, initial state, constraints,
and battery use are comparable. Constraint violations and solver failures remain
visible in the result tables.
