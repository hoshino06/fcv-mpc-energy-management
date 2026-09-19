function S = build_udds_scenario()
%BUILD_UDDS_SCENARIO Single source of truth for the UDDS study.
%
% S.udds_file, S.dt, S.traction_scale : demand construction
% S.plant, S.limits, S.solver         : shared plant / constraints / proposed solver
% S.mt_min_stop_s, S.mt_pad_s         : micro-trip segmentation parameters
% S.select                            : which micro-trips to run
%                                       'all' | vector of indices | 'representative'
% S.segments(k)                       : one micro-trip (from micro_trips.m) plus
%                                       .name .Qmax_As .lowpass_tau_s
%
% Qmax_As is left [] on purpose: run_study sets it per segment to the net
% battery discharge the frequency-decoupling method actually used, so both
% methods are compared at equal battery utilisation (paper Sec. V-B).

cfgdir = fileparts(mfilename('fullpath'));                 % .../experiments/studies/03_udds
expdir = fileparts(fileparts(cfgdir));                     % .../experiments
S.study = '03_udds';
S.udds_file = fullfile(cfgdir,'data', ...
                       '都市ダイナモメーター運転スケジュール_PFt付き.txt');
if ~isfile(S.udds_file)
    a = dir(fullfile(expdir,'**','都市ダイナモメーター*PFt*.txt'));
    if ~isempty(a), S.udds_file = fullfile(a(1).folder,a(1).name); end
end

S.dt             = 0.05;
S.plant_dt       = 0.005; % common RK4 plant step, independent of controller period
S.traction_scale = 1/5;

% --- shared plant / initial condition ----------------------------------
S.plant.start_Vcm_idx = 100;
S.plant.start_Ifc_idx = 100;
S.plant.start_SOC     = 0.5;
S.plant.battery_parallel = 15;
S.plant.battery_series   = 12;

% --- constraint limits (post-hoc, identical for every method) ----------
% All logged currents/voltages are PHYSICAL: Vcm [V], Ist/Ib [A].
S.limits.lambda_O2_min = 1.5;
S.limits.Ist_min       = 0;      % A
S.limits.Ist_max       = 616;    % A  (2.2 * 280)
S.limits.Vcm_min       = 0;      % V
S.limits.Vcm_max       = 300;    % V  (params.Vcm_u)
S.limits.Ib_abs_max    = 36;     % A
S.limits.Qmax_tol_As   = 0.5;

% --- proposed-method headline solver ---------------------------------
S.solver.use_rk4      = 1;
S.solver.integ_h      = 0.025;
S.solver.ls_backtrack = 1;
S.solver.N_iter       = 50;
S.solver.N_iter_ref   = 200;

% --- segmentation: micro-trips (context) then fixed T-second intervals ---
S.mt_min_stop_s  = 5;     % split UDDS on stops >= this long
S.mt_pad_s       = 2;     % micro-trip settling pad
S.mt_max_trip_s  = 90;    % cut micro-trips longer than this at interior speed minima
S.interval_T_s        = 5;    % IDC interval length: fixed tiling of each micro-trip's moving span
S.interval_pad_s      = 1;    % simulate this much past the scored window
S.eval_settle_s       = 0;    % discard this much at the start in eval. 0 = score the
                              % whole window (on short from-equilibrium intervals the
                              % initial step response IS the event and must be kept).
S.interval_min_peak_kW = 3;   % intervals below this scaled peak power are "idle" -> not an IDC interval

% which intervals to run:  'representative' | 'all' (= active) | [i j ...] | 'mtNN'
S.select = 'representative';

mt = micro_trips(S.udds_file, S.mt_min_stop_s, S.mt_pad_s, S.mt_max_trip_s);
iv = fixed_intervals(S.udds_file, mt, S.interval_T_s, S.interval_pad_s);

S.all_micro_trips = mt;
S.all_intervals   = iv;
active = ([iv.peak_kW] * S.traction_scale) >= S.interval_min_peak_kW;
S.active_idx = find(active);           % indices into S.all_intervals that carry real load

% representative = terciles of transient severity (scaled peak power) over ACTIVE intervals
act = S.active_idx;
sev = [iv(act).peak_kW] * S.traction_scale;
[~,ord] = sort(sev);
q = @(f) act(ord(max(1,round(f*numel(ord)))));
S.rep_idx = unique([q(1/6) q(1/2) q(5/6)]);

if ischar(S.select)
    if strcmp(S.select,'all')
        sel = S.active_idx;
    elseif startsWith(S.select,'mt')
        want = sscanf(S.select(3:end),'%d');
        sel = intersect(find([iv.mt_idx] == want), S.active_idx);
    else
        sel = S.rep_idx;
    end
else
    sel = S.select;
end

S.segments = struct('name',{},'t0',{},'t1',{},'sim_t1',{},'peak_kW',{}, ...
                    'energy_kJ',{},'Qmax_As',{},'lowpass_tau_s',{}, ...
                    'mt_idx',{},'mt_part',{},'iv_idx',{});
for j = 1:numel(sel)
    g = iv(sel(j));
    S.segments(j) = struct( ...
        'name', sprintf('iv%03d_mt%d%s', g.idx, g.mt_idx, g.mt_part), ...
        't0', g.t0, 't1', g.t1, 'sim_t1', g.sim_t1, ...
        'peak_kW', g.peak_kW, 'energy_kJ', g.energy_kJ, ...
        'Qmax_As', [], 'lowpass_tau_s', [], ...
        'mt_idx', g.mt_idx, 'mt_part', g.mt_part, 'iv_idx', g.idx);
end
end
