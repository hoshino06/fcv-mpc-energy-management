function iv = fixed_intervals(udds_file, mt, T_s, pad_s)
%FIXED_INTERVALS  Tile the moving span of each micro-trip with fixed T-second intervals.
%
%   iv = fixed_intervals(udds_file, mt, T_s, pad_s)
%
%   mt     : micro-trip struct array from micro_trips.m
%   T_s    : interval length [s] (default 5)
%   pad_s  : seconds of settle kept AFTER each interval's window for the sim
%            (the plant/optimizer needs a little past the scored span); default 1
%
% Each micro-trip's moving portion [move_t0, move_t1] is divided into
% consecutive T_s windows; a trailing remainder shorter than T_s/2 is merged
% into the previous window. Every interval is an independent IDC interval with
% its own discharge budget.
%
% Returns iv(k) with:
%   .idx        global interval number
%   .mt_idx .mt_part
%   .t0 .t1     scored window [s]
%   .sim_t1     t1 + pad_s (simulate to here, score to t1)
%   .peak_kW .energy_kJ  from the unscaled PFt column over [t0,t1]
%   .in_mt_pos  1..N position within its micro-trip

if nargin < 3 || isempty(T_s),   T_s = 5;   end
if nargin < 4 || isempty(pad_s), pad_s = 1; end

T = readtable(udds_file,'VariableNamingRule','preserve');
tt  = T{:,1};
PFt = T{:,end};
dt  = median(diff(tt));

iv = struct('idx',{},'mt_idx',{},'mt_part',{},'t0',{},'t1',{},'sim_t1',{}, ...
            'peak_kW',{},'energy_kJ',{},'in_mt_pos',{});
g = 0;
for m = 1:numel(mt)
    a = mt(m).move_t0;  b = mt(m).move_t1;
    edges = a:T_s:b;
    if numel(edges) < 2, edges = [a b]; end
    if b - edges(end) >= T_s/2, edges(end+1) = b; else, edges(end) = b; end %#ok<AGROW>
    for p = 1:numel(edges)-1
        w0 = edges(p);  w1 = edges(p+1);
        sel = tt >= w0 & tt <= w1;
        g = g + 1;
        s = struct();
        s.idx = g;
        s.mt_idx = mt(m).idx;  s.mt_part = mt(m).part;
        s.t0 = w0;  s.t1 = w1;
        s.sim_t1 = min(tt(end), w1 + pad_s);
        s.peak_kW   = max(PFt(sel))/1e3;
        s.energy_kJ = sum(PFt(sel))*dt/1e3;
        s.in_mt_pos = p;
        iv(end+1) = s; %#ok<AGROW>
    end
end
end
