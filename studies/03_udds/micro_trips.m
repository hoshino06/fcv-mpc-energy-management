function mt = micro_trips(udds_file, min_stop_s, pad_s, max_trip_s)
%MICRO_TRIPS  Segment a drive cycle into micro-trips (stop-to-stop).
%
%   mt = micro_trips(udds_file, min_stop_s, pad_s, max_trip_s)
%
% A micro-trip is the interval between two sustained stops. A "sustained stop"
% is a run of >= min_stop_s consecutive samples at zero speed. This is the
% standard drive-cycle decomposition used in EMS / drive-cycle analysis and
% removes any subjective window selection.
%
%   min_stop_s : minimum stop duration to split on (default 5 s)
%   pad_s      : seconds of the leading stop kept before motion, and of the
%                trailing stop kept after (settling); default 2 s
%   max_trip_s : if > 0, micro-trips whose moving span exceeds this are cut at
%                interior local speed minima until every piece is <= max_trip_s
%                (default 0 = no splitting). Split pieces share .idx and get a
%                .part letter.
%
% Returns a struct array mt(k) with fields:
%   .idx   micro-trip number      .part  '' or 'a','b',...
%   .t0 .t1     UDDS time window [s] (integer, includes pad)
%   .dur_s      t1 - t0
%   .move_t0 .move_t1   the moving portion (no pad)
%   .peak_kW .energy_kJ  from the (unscaled) PFt column
%   .n_stops_skipped     short stops (< min_stop_s) inside this micro-trip

if nargin < 2 || isempty(min_stop_s), min_stop_s = 5; end
if nargin < 3 || isempty(pad_s),      pad_s = 2;      end
if nargin < 4 || isempty(max_trip_s), max_trip_s = 0; end

T = readtable(udds_file,'VariableNamingRule','preserve');
t  = T{:,1};
% speed column: prefer m/s (col 3), fall back to mph (col 2)
if size(T,2) >= 3, v = T{:,3}; else, v = T{:,2}; end
PFt = T{:,end};
dt  = median(diff(t));                       % 1 s for UDDS

moving = v > 0.05;                            % m/s threshold
% runs of not-moving
nm = ~moving;
d  = diff([0; nm; 0]);
stop_start = find(d==1);
stop_end   = find(d==-1) - 1;
stop_dur   = (stop_end - stop_start + 1) * dt;
long_stop  = stop_dur >= min_stop_s;
ss = stop_start(long_stop);  se = stop_end(long_stop);

% micro-trip k spans from end of long-stop k to start of long-stop k+1
mt = struct('idx',{},'part',{},'t0',{},'t1',{},'dur_s',{},'move_t0',{},'move_t1',{}, ...
            'peak_kW',{},'energy_kJ',{},'n_stops_skipped',{});
np = round(pad_s/dt);
mtnum = 0;
for k = 1:numel(ss)-1
    a_move = se(k) + 1;                       % first moving sample
    b_move = ss(k+1) - 1;                     % last moving sample
    if b_move <= a_move, continue; end
    mtnum = mtnum + 1;
    short_inside = ~long_stop & stop_start > a_move & stop_end < b_move;
    nskip = nnz(short_inside);

    % split the moving span [a_move,b_move] at interior local speed minima
    cuts = split_span(a_move, b_move, v, max_trip_s, dt);
    nparts = numel(cuts) - 1;
    for p = 1:nparts
        am = cuts(p); bm = cuts(p+1);
        if p > 1,      am = am + 1; end       % avoid overlap
        a = max(1, am - np*(p==1));
        b = min(numel(t), bm + np*(p==nparts));
        seg = PFt(am:bm);
        g = struct();
        g.idx  = mtnum;
        g.part = ''; if nparts > 1, g.part = char('a'+p-1); end
        g.t0 = t(a);           g.t1 = t(b);
        g.dur_s = g.t1 - g.t0;
        g.move_t0 = t(am);     g.move_t1 = t(bm);
        g.peak_kW   = max(seg)/1e3;
        g.energy_kJ = sum(seg)*dt/1e3;
        g.n_stops_skipped = nskip;
        mt(end+1) = g; %#ok<AGROW>
    end
end
end

function cuts = split_span(a, b, v, max_trip_s, dt)
% recursively cut [a,b] (sample indices) at the interior local speed minimum
% until each piece's duration <= max_trip_s. Returns sorted sample-index cuts.
cuts = [a b];
if max_trip_s <= 0, return; end
changed = true;
while changed
    changed = false;
    newcuts = cuts(1);
    for i = 1:numel(cuts)-1
        p = cuts(i); q = cuts(i+1);
        if (q - p)*dt <= max_trip_s
            newcuts(end+1) = q; %#ok<AGROW>
            continue;
        end
        % candidate cut = index of minimum speed in the interior third..two-thirds
        lo = p + round((q-p)*0.3);  hi = p + round((q-p)*0.7);
        [~,rel] = min(v(lo:hi));
        c = lo + rel - 1;
        newcuts(end+1) = c; %#ok<AGROW>
        newcuts(end+1) = q; %#ok<AGROW>
        changed = true;
    end
    cuts = unique(newcuts);
end
end
