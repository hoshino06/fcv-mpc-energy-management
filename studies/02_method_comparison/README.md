# 02 — Conventional energy-management comparison

This second paper stage compares the proposed DDP-based controller with a
low-pass frequency-decoupling benchmark and a direct-shooting fmincon reference
on a common two-step demand. The demand increments are reduced to 8 kW to avoid
the steady-state boundary encountered by the 10 kW concept demonstration.

```matlab
addpath('simulation/fcv-mpc-energy-management/studies/02_method_comparison');
run_study
```

The runner first evaluates the low-pass benchmark at `omega_n=16`, then sets
the NMPC discharge budget to the measured low-pass battery use. All methods use
the common 12-series by 45-parallel equivalent battery and a +/-36 A pack-current
limit where that constraint is implemented.

Results are written to `results/02_method_comparison/`. The current data are
diagnostic, not a final paper comparison: the low-pass benchmark does not enforce
the common oxygen-excess-ratio or battery-current limits, and the nonuniform N=5
case is not yet converged to the matched discharge budget. These limitations are
retained in `method_comparison_8kw.csv`.
