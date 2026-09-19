# Step-demand experiment

This scenario recovers the two-step demand used in Section V-A (Figs. 3--6)
of `FCV_MPC.pdf`: the equilibrium demand increases by 10 kW at 2.0 s and by
another 10 kW at 2.2 s. Iwai's setup constructs 20 s of preview data, while
the simulations and hydrogen totals used in Figs. 3--6 cover 0--5 s.

The numerical model and solver are shared with `experiments/methods/ddp`.
Only the demand and experiment parameters live here; execution has no runtime
dependency on `ronbun_iwai`.

For a less boundary-limited comparison, use `step_demand_scenario(8000)`;
this applies two 8 kW increments and gives a final demand of about 39.42 kW.
The comparison uses a common equivalent battery of 12 series x 45 parallel
cells with a +/-36 A pack-current limit. This is a numerically robust test-bed
configuration, not a claimed production-FCV cell topology.
`run_method_comparison` first runs the low-pass benchmark at `omega_n=16`, then
sets the NMPC discharge budget to its measured battery use. The battery dynamics
are regenerated for the common 12-series x 45-parallel equivalent pack. Pack
current and per-cell current must still be distinguished in reporting.

The current matched-use data are diagnostic rather than the final paper result.
The low-pass benchmark tracks power well but does not enforce the common oxygen
excess-ratio or +/-36 A NMPC bounds. The nonuniform N=5 run with 20 iterations is
real-time but is not converged to the discharge budget; 50 iterations drives the
plant model outside the compressor map. These limitations are written into the
comparison CSV rather than hidden.

```matlab
addpath('simulation/experiments/studies/step_demand');
run_study
```

The study uses `dt=0.05 s`, uniform `N=10` (`Tpred=0.5 s`), 100 DDP
iterations, and `Qmax={72,36,18,3.6} As`. The index, summary, and figures are
written to `simulation/experiments/results/step_demand`; raw runs are stored in
`results/step_demand/ddp`.
