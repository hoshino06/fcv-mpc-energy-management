function M = eval_run(resfile, limits_override)
%EVAL_RUN  Method-agnostic metrics from one saved benchmark run.
%
%   M = eval_run('.../results/<method>/<tag>.mat')
%
% The .mat must contain: simout, meta, params, Ref, seg, limits, demand.
%
% IMPORTANT: logged signals come out at the Simulink solver's own timestep,
% which is NOT params.MPC.dt and differs between runs. All windowing here uses
% each signal's real .Time vector. Accumulators (H2, S) are read by
% interpolating to the window edges.

R  = load(resfile);
so = R.simout; params = R.params; Ref = R.Ref; seg = R.seg; L = R.limits;
M  = struct(); M.notes = {}; M.tag = R.meta.tag; M.method = R.meta.method;
if nargin>=2 && ~isempty(limits_override)
    L=limits_override;
    M.notes{end+1}='constraints evaluated with explicit common physical limits';
elseif L.Vcm_max==1 && Ref.Vcm==300
    M.notes{end+1}='legacy Vcm bound is normalized; pass physical limits as second argument';
end

sw = [0, inf];
if isfield(R.meta,'score_win_s') && numel(R.meta.score_win_s)==2, sw = R.meta.score_win_s(:).'; end

    function v = getsig(name)              % raw timeseries/array from simout or a struct
        v = [];
        if isstruct(so)
            if isfield(so,name), v = so.(name); end
        else
            try, v = so.get(name); catch, end
        end
    end

    function [d,tk] = sig(name)             % data + time, trimmed to the window
        d = []; tk = [];
        v = getsig(name);
        if isempty(v), return; end
        if isa(v,'timeseries')
            d = v.Data; tk = v.Time(:);
        elseif isnumeric(v)
            d = v; tk = (0:size(v,1)-1).' * params.MPC.dt;
        else, return; end
        d = squeeze(d);
        if isvector(d), d = d(:); end
        if size(d,1) ~= numel(tk) && size(d,2) == numel(tk), d = d.'; end
        keep = tk >= sw(1) - 1e-9 & tk <= sw(2) + 1e-9;
        d = d(keep,:); tk = tk(keep);
    end
    function y = acc_delta(name)            % accumulator increase across the window
        y = [];
        v = getsig(name);
        if ~isa(v,'timeseries'), return; end
        dd = squeeze(v.Data); dd = dd(:); tt = v.Time(:);
        lo = interp1(tt, dd, max(sw(1), tt(1)),  'linear', 'extrap');
        hi = interp1(tt, dd, min(sw(2), tt(end)),'linear', 'extrap');
        y = hi - lo;
    end
    function y = acc_max(name)              % max increase from window start
        y = [];
        v = getsig(name);
        if ~isa(v,'timeseries'), return; end
        dd = squeeze(v.Data); dd = dd(:); tt = v.Time(:);
        m  = tt >= sw(1) - 1e-9 & tt <= sw(2) + 1e-9;
        if ~any(m), return; end
        lo = interp1(tt, dd, max(sw(1), tt(1)), 'linear', 'extrap');
        y  = max(dd(m)) - lo;
    end

% ---------- timing ----------
tI = sig('time_ILQR');
if isempty(tI)
    M.notes{end+1} = 'no time_ILQR (no per-step timing for this method yet)';
else
    v = tI(isfinite(tI) & tI > 0);
    if numel(v) > 1, v = v(2:end); end
    M.n_updates   = numel(v);
    M.t_median_ms = 1e3*median(v);
    M.t_p95_ms    = 1e3*prctile_(v,95);
    M.t_max_ms    = 1e3*max(v);
    M.deadline_ms = 1e3*params.MPC.dt;
    M.miss_rate   = mean(v > params.MPC.dt);
end

% ---------- oxygen excess ratio ----------
z = sig('z'); if isempty(z), z = sig('sannso'); end
if ~isempty(z)
    z = z(:); z = z(isfinite(z));
    M.lambdaO2_min  = min(z);
    M.lambdaO2_viol = max(0, L.lambda_O2_min - M.lambdaO2_min);
else
    M.notes{end+1} = 'no z / sannso signal';
end

% ---------- currents (all models log Vcm/Ist/Ib in PHYSICAL units: V, A, A) ----------
Ifc = sig('Ifc');                                 % stack current [A], bound [0, 616]
if ~isempty(Ifc)
    M.Ist_max = max(Ifc(:)); M.Ist_min = min(Ifc(:));
    M.Ist_viol = max([0; Ifc(:) - L.Ist_max; L.Ist_min - Ifc(:)]);
end
Ib = sig('Ib');                                   % battery current [A], bound +/- 36
if ~isempty(Ib)
    M.Ib_abs_max = max(abs(Ib(:)));
    M.Ib_viol    = max(0, M.Ib_abs_max - L.Ib_abs_max);
end
Vcm = sig('Vcm');                                 % compressor voltage [V], bound [0, 300]
if ~isempty(Vcm)
    M.Vcm_max = max(Vcm(:)); M.Vcm_min = min(Vcm(:));
    M.Vcm_viol = max([0; Vcm(:) - L.Vcm_max; L.Vcm_min - Vcm(:)]);
end

% ---------- battery use over the scored window ----------
dS = acc_delta('S');
if ~isempty(dS)
    M.batt_used_As     = dS * Ref.Ib;
    M.batt_used_As_max = acc_max('S') * Ref.Ib;
    if isfield(seg,'Qmax_As') && ~isempty(seg.Qmax_As)
        M.Qmax_As           = seg.Qmax_As;
        M.Qmax_overshoot_As = max(0, M.batt_used_As_max - seg.Qmax_As);
        M.Qmax_use_frac     = M.batt_used_As / seg.Qmax_As;
    end
else
    M.notes{end+1} = 'no S accumulator';
end

% ---------- SOC drift over the window ----------
soc = sig('Battery_state'); if isempty(soc), soc = sig('SOC'); end
if ~isempty(soc)
    s1 = soc(1,1); s2 = soc(end,1);
    M.SOC_start = s1; M.SOC_end = s2; M.SOC_drop = s1 - s2;
    M.batt_from_SOC_As = M.SOC_drop * 3060 * params.battery.parallel;
    if isfield(M,'batt_used_As')
        M.batt_SOC_discrepancy_As=M.batt_from_SOC_As-M.batt_used_As;
    end
end

% ---------- hydrogen over the scored window ----------
dH2 = acc_delta('H2');
if ~isempty(dH2) && isfield(params,'n')
    M.H2_g = dH2 * 2.016 * params.n * Ref.Ifc / (2*96485.3);
elseif ~isempty(dH2)
    M.H2_raw_win = dH2;
    M.notes{end+1} = 'params.n missing -> H2 raw';
else
    M.notes{end+1} = 'no H2 accumulator';
end

% ---------- power tracking: NET delivered power vs demand ----------
% "net delivered" = Pfc - Pcm + Pb (what the OCP tracks). Every model logs it
% as `Pd`; only fall back to reconstructing it.
[deliv,tk] = net_delivered_power(R);
keep = tk>=sw(1)-1e-9 & tk<=sw(2)+1e-9;
deliv=deliv(keep); tk=tk(keep);
if ~isempty(deliv) && isfield(R,'demand') && isa(R.demand,'timeseries')
    deliv = deliv(:);
    dref  = interp1(R.demand.Time, R.demand.Data(:,1), tk, 'linear', 'extrap');
    e = deliv - dref;
    M.track_RMSE_W   = sqrt(trapz(tk,e.^2)/(tk(end)-tk(1)));
    M.track_maxabs_W = max(abs(e));
    M.deliv_mean_W   = trapz(tk,deliv)/(tk(end)-tk(1));
    M.demand_mean_W  = trapz(tk,dref)/(tk(end)-tk(1));
else
    M.notes{end+1} = 'power-tracking not computed';
end
ef=sig('exitflag');
if ~isempty(ef), M.exitflag_min=min(ef); M.solver_failed_updates=sum(ef<=0); end
fo=sig('firstorderopt'); if ~isempty(fo), M.firstorderopt_max=max(fo); end
end

function q = prctile_(v,p)
v = sort(v(:)); n = numel(v);
if n==0, q=NaN; return; end
if n < 2, q = v(end); return; end
q = interp1((0.5:1:n-0.5)/n*100, v, p, 'linear', 'extrap');
q = min(max(q, v(1)), v(end));
end
