function figfile = plot_interval(seg_name, tags)
%PLOT_INTERVAL  Overlay the time-domain waveforms of several methods on one interval.
%
%   plot_interval('iv208_mt13')                 % auto: every results/*/*iv208_mt13*.mat
%   plot_interval('iv208_mt13', {'ddp/iv208_mt13_nonunif5_N50', ...
%                                'ddp/iv208_mt13_unif10_N50', ...
%                                'fmincon/iv208_mt13_unifNp10', ...
%                                'lowpath/iv208_mt13_wn16'})
%
% Panels: delivered vs demanded power, Pfc, Pb, lambda_O2 (+1.5), SOC/qdis,
% Ist(+bounds), Ib(+/-1 norm), Vcm, MPC solve time (+50 ms).
% -> results/waveforms/<seg>.png

studydir = fileparts(mfilename('fullpath'));
expdir   = fileparts(fileparts(studydir));
resroot  = fullfile(expdir,'results','03_udds');

if nargin < 2 || isempty(tags)
    f = dir(fullfile(resroot,'*',['*' seg_name '*.mat']));
    files = arrayfun(@(x) fullfile(x.folder,x.name), f, 'uni', 0);
else
    files = cell(size(tags));
    for i = 1:numel(tags)
        t = tags{i};
        if isfile(t), files{i} = t;
        elseif isfile(fullfile(resroot,[t '.mat'])), files{i} = fullfile(resroot,[t '.mat']);
        else, files{i} = fullfile(resroot,t); end
    end
end
files = files(cellfun(@isfile, files));
assert(~isempty(files), 'no result files for "%s"', seg_name);

R1 = load(files{1});
sw = [0 inf]; if isfield(R1.meta,'score_win_s'), sw = R1.meta.score_win_s(:).'; end

% panel = {signal-key, ylabel, scale, hlines, ylog}
P = { 'DELIV', 'P deliv/dem [kW]', 1e-3, [],       false;
      'Pb',    'P_{b} [kW]',    1e-3, 0,        false;
      'Z',     '\lambda_{O2}',  1,    1.5,      false;
      'SOC',   'SOC',           1,    [],       false;
      'S',     'q_{dis} [norm]',1,    [],       false;
      'Ifc',   'I_{st} [A]',    1,    [0 616],  false;
      'Ib',    'I_b [A]',       1,    [-36 36], false;
      'Vcm',   'V_{cm} [V]',    1,    [0 300],  false;
      'time_ILQR','solve [ms]', 1e3,  50,       true };
np = size(P,1);

fig = figure('Name',['waveforms ' seg_name],'Position',[40 30 1150 1600]);
tl  = tiledlayout(np,1,'TileSpacing','compact');
ax  = gobjects(np,1);
for p = 1:np, ax(p) = nexttile; hold(ax(p),'on'); grid(ax(p),'on'); ylabel(ax(p),P{p,2}); end
cols = lines(max(7,numel(files)));
leg  = strings(1,numel(files));
hmet = gobjects(1,numel(files));       % one line handle per method (from panel 1)

hdem = plot(ax(1), R1.demand.Time, R1.demand.Data(:,1)*1e-3, 'k--', 'LineWidth', 1.8);

for i = 1:numel(files)
    R = load(files{i});
    [pdir,fn] = fileparts(fileparts(files{i}));  %#ok<ASGLU>
    [~,mth]   = fileparts(fileparts(files{i}));
    [~,base]  = fileparts(files{i});
    leg(i) = string(mth) + "/" + extractAfter(string(base), string(seg_name)+"_");
    c = cols(i,:);
    for p = 1:np
        key = P{p,1};
        switch key
            case 'DELIV'
                [d,td] = net_delivered_power(R);
            case 'Z'
                [d,td] = getsig(R,'z'); if isempty(d), [d,td] = getsig(R,'sannso'); end
            case 'SOC'
                [d,td] = getsig(R,'Battery_state'); if isempty(d), [d,td] = getsig(R,'SOC'); end
                if ~isempty(d) && size(d,2) > 1, d = d(:,1); end
            otherwise
                [d,td] = getsig(R,key);
        end
        if isempty(d), continue; end
        h = plot(ax(p), td(:), d(:)*P{p,3}, '-', 'Color', c, 'LineWidth', 1.1);
        if p == 1, hmet(i) = h; end
    end
end

% limit lines
for p = 1:np
    L = P{p,4};
    for v = L(:).', yline(ax(p), v, 'r:', 'LineWidth', 1); end
    xline(ax(p), sw(2), 'Color',[.6 .6 .6]);
    if P{p,5}, set(ax(p),'YScale','log'); end
end
xlabel(ax(np),'time in interval [s]');
hh = [hdem, hmet(isgraphics(hmet))];
ll = ["demand", leg(isgraphics(hmet))];
legend(ax(1), hh, ll, 'Interpreter','none', 'Location','eastoutside');
title(tl, sprintf('%s   (scored [%.2g, %.2g] s)', seg_name, sw(1), sw(2)), 'Interpreter','none');

outdir = fullfile(resroot,'waveforms'); if ~isfolder(outdir), mkdir(outdir); end
figfile = fullfile(outdir,[seg_name '.png']);
saveas(fig, figfile);
fprintf('wrote %s  (%d methods)\n', figfile, numel(files));
end

% -------------------------------------------------------------------------
function [d,tk] = getsig(R, name)
d = []; tk = [];
so = R.simout;
if isstruct(so)
    if isfield(so,name), v = so.(name); else, return; end
else
    try, v = so.get(name); catch, return; end
end
if isa(v,'timeseries')
    d = squeeze(v.Data); tk = v.Time(:);
elseif isnumeric(v)
    d = squeeze(v); tk = (0:size(v,1)-1).' * R.params.MPC.dt;
else
    return;
end
if isvector(d), d = d(:); end
if size(d,1) ~= numel(tk) && size(d,2) == numel(tk), d = d.'; end
end
