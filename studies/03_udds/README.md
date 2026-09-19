# 03 — UDDS validation

This study compares DDP with uniform/non-uniform horizons, fmincon, and the
low-pass frequency-decoupling method on a predeclared UDDS interval.

```matlab
addpath('simulation/experiments/studies/03_udds');
run_study                 % all phases
run_study('lowpath')      % frequency-decoupling runs
run_study('nmpc')         % DDP and fmincon runs
run_study('report')       % table and waveform figure
```

`udds_scenario.m` is the public configuration entry point.
The source demand file is kept in `data/` beside the UDDS-specific scripts.
Results are written to `results/03_udds/`. This is the third paper stage; its
final manuscript framing is still open.
