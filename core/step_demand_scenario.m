function S = step_demand_scenario(step_increment_W, study_name)
%STEP_DEMAND_SCENARIO Shared step-demand configuration for studies 01 and 02.
if nargin<1, step_increment_W=10000; end
if nargin<2, study_name='01_qmax_sweep'; end
S.study=study_name;
S.dt=0.05;
S.plant_dt=0.005;
S.traction_scale=1/5;
S.plant=struct('start_Vcm_idx',100,'start_Ifc_idx',100,'start_SOC',0.5, ...
    'battery_parallel',45,'battery_series',12);
S.limits=struct('lambda_O2_min',1.5,'Ist_min',0,'Ist_max',616, ...
    'Vcm_min',0,'Vcm_max',300,'Ib_abs_max',36,'Qmax_tol_As',0.5);
S.solver=struct('use_rk4',1,'integ_h',0.025,'ls_backtrack',1, ...
    'N_iter',50,'N_iter_ref',200);
S.eval_settle_s=0;
S.demand.type='step';
S.demand.step_times_s=[2.0 2.2];
S.demand.step_increments_W=[step_increment_W step_increment_W];
S.Tpred_s=0.5;
S.figure_time_s=5;
S.Qmax_As=[72 36 18 3.6];
S.segments=struct('name','step_demand','t0',0,'t1',5,'sim_t1',5, ...
    'peak_kW',NaN,'energy_kJ',NaN,'Qmax_As',72,'lowpass_tau_s',[], ...
    'mt_idx',NaN,'mt_part','','iv_idx',NaN);
end
