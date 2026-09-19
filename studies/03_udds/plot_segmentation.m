%% plot_segmentation.m
% 論文の「区間の選び方」説明図（3 面）:
%   (1) UDDS 速度全体 + マイクロトリップ網掛け（外側構造）
%   (2) 1 マイクロトリップを拡大：電力トレース + 5 秒固定区間グリッド + 区間番号
%   (3) 全 5 秒区間の重大度（スケール済みピーク電力）ヒストグラム + 代表区間マーク
%
% UDDS study directoryをMATLAB pathに追加して実行する。

S = build_udds_scenario();
T = readtable(S.udds_file,'VariableNamingRule','preserve');
t   = T{:,1};
v   = T{:,3};
PFt = T{:,end} * S.traction_scale / 1e3;             % scaled traction power [kW]
mt  = S.all_micro_trips;
iv  = S.all_intervals;
act = S.active_idx;
rep = S.rep_idx;
ivpk = [iv.peak_kW] * S.traction_scale;

fprintf('\nUDDS: %d micro-trips (stop>=%gs, cut>%gs) -> %d fixed %g s intervals, %d active (peak >= %g kW)\n', ...
    numel(mt), S.mt_min_stop_s, S.mt_max_trip_s, numel(iv), S.interval_T_s, ...
    numel(act), S.interval_min_peak_kW);
moving_s = sum(arrayfun(@(g) g.t1-g.t0, iv(act)));
fprintf('scored active driving time = %.0f s of %.0f s cycle (%.0f%%)\n', ...
    moving_s, t(end), 100*moving_s/t(end));
fprintf('representative intervals (peak-power terciles):\n');
for r = rep
    g = iv(r);
    fprintf('  iv%03d  mt%d%s pos %d  t %d-%d  peak %.1f kW  E %.1f kJ\n', ...
        g.idx, g.mt_idx, g.mt_part, g.in_mt_pos, g.t0, g.t1, ...
        g.peak_kW*S.traction_scale, g.energy_kJ*S.traction_scale);
end

fig = figure('Name','UDDS segmentation','Position',[60 60 1400 780]);
tiledlayout(3,1,'TileSpacing','compact');

% (1) full speed + micro-trips
ax1 = nexttile; hold(ax1,'on'); grid(ax1,'on');
plot(ax1,t,v,'Color',[0.35 0.35 0.35]);
ylabel(ax1,'speed [m/s]'); yl = ylim(ax1);
for k = 1:numel(mt)
    g = mt(k);
    patch(ax1,[g.t0 g.t1 g.t1 g.t0],[yl(1) yl(1) yl(2) yl(2)], [0.2 0.5 1], ...
        'FaceAlpha',0.10,'EdgeColor','none');
    text(ax1,mean([g.t0 g.t1]),yl(2)*0.9,sprintf('%d%s',g.idx,g.part), ...
        'HorizontalAlignment','center','FontSize',7);
end
title(ax1,sprintf('(1) UDDS split into %d micro-trips at stops \\geq %g s', numel(mt), S.mt_min_stop_s));
xlim(ax1,[0 t(end)]);

% (2) zoom on the micro-trip that holds the middle representative
gz = iv(rep(ceil(numel(rep)/2)));
mzi = find([mt.idx]==gz.mt_idx & strcmp({mt.part},gz.mt_part),1);
mz  = mt(mzi);
ax2 = nexttile; hold(ax2,'on'); grid(ax2,'on');
sel = t>=mz.t0 & t<=mz.t1;
plot(ax2,t(sel),PFt(sel),'k-','LineWidth',1.2);
ylabel(ax2,'scaled traction power [kW]'); yl2 = ylim(ax2);
ivm = iv([iv.mt_idx]==mz.idx & strcmp({iv.mt_part},mz.part));
for j = 1:numel(ivm)
    g = ivm(j);
    xline(ax2,g.t0,'-','Color',[0.6 0.6 0.6]);
    isrep = any(rep==g.idx);
    if isrep
        patch(ax2,[g.t0 g.t1 g.t1 g.t0],[yl2(1) yl2(1) yl2(2) yl2(2)],[1 0.6 0.1], ...
            'FaceAlpha',0.25,'EdgeColor','none');
    end
    fw = 'normal'; if isrep, fw = 'bold'; end
    text(ax2,mean([g.t0 g.t1]),yl2(2)*0.92,sprintf('%d',g.in_mt_pos), ...
        'HorizontalAlignment','center','FontSize',8,'FontWeight',fw);
end
xline(ax2,mz.t1,'-','Color',[0.6 0.6 0.6]);
title(ax2,sprintf('(2) micro-trip %d%s tiled into %g s IDC intervals (orange = a representative)', ...
    mz.idx,mz.part,S.interval_T_s));
xlabel(ax2,'UDDS time [s]'); xlim(ax2,[mz.t0 mz.t1]);

% (3) severity histogram over ACTIVE intervals
ax3 = nexttile; hold(ax3,'on'); grid(ax3,'on');
histogram(ax3,ivpk(act),20,'FaceColor',[0.4 0.6 0.9]);
for r = rep
    xline(ax3,ivpk(r),'r-','LineWidth',1.5, ...
        'Label',sprintf('iv%03d',iv(r).idx),'LabelOrientation','horizontal');
end
xlabel(ax3,'interval peak scaled traction power [kW]');
ylabel(ax3,'# intervals');
title(ax3,sprintf('(3) severity of %d active intervals; representatives = peak-power terciles', numel(act)));

studydir=fileparts(mfilename('fullpath'));
outdir = fullfile(fileparts(fileparts(studydir)),'results','03_udds');
if ~isfolder(outdir), mkdir(outdir); end
saveas(fig, fullfile(outdir,'udds_segmentation.png'));
fprintf('\nfigure -> %s\n', fullfile(outdir,'udds_segmentation.png'));
