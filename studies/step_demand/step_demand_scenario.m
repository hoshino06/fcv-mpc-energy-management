function S = step_demand_scenario(step_increment_W)
%STEP_DEMAND_SCENARIO Paper Sec. V-A demand in the common benchmark format.
if nargin<1, step_increment_W=10000; end
here=fileparts(mfilename('fullpath')); expdir=fileparts(fileparts(here));
addpath(fullfile(expdir,'studies','udds'));
S=build_udds_scenario();
S.study='step_demand';
S.demand.type='step';
S.demand.step_times_s=[2.0 2.2];
S.demand.step_increments_W=[step_increment_W step_increment_W];
% Numerically robust common equivalent pack for the method comparison.
S.plant.battery_parallel=45;
S.plant.battery_series=12;
S.limits.Ib_abs_max=36;
S.Tpred_s=0.5;
S.figure_time_s=5;
S.Qmax_As=[72 36 18 3.6];
S.segments=struct('name','step_demand','t0',0,'t1',5,'sim_t1',5, ...
    'peak_kW',NaN,'energy_kJ',NaN,'Qmax_As',72,'lowpass_tau_s',[], ...
    'mt_idx',NaN,'mt_part','','iv_idx',NaN);
end
