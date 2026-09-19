%% collect_table.m
% results/03_udds/<method>/*.mat を評価し、UDDSの統合比較表を作る。

here   = fileparts(mfilename('fullpath'));
expdir = fileparts(fileparts(here));
addpath(here,fullfile(expdir,'core'));
common_scenario=build_udds_scenario();
resroot = fullfile(expdir,'results','03_udds');
files=[dir(fullfile(resroot,'ddp','*.mat'));dir(fullfile(resroot,'fmincon','*.mat'));dir(fullfile(resroot,'lowpath','*.mat'))];
if isempty(files), fprintf('no result files under %s\n', resroot); return; end

rows = struct([]);
for i = 1:numel(files)
    f = fullfile(files(i).folder, files(i).name);
    try
        M = eval_run(f,common_scenario.limits);
    catch ME
        fprintf('  eval failed for %s: %s\n', files(i).name, ME.message);
        continue;
    end
    g = struct();
    [~, seg_from_dir] = fileparts(files(i).folder);
    g.method   = string(getf(M,'method',seg_from_dir));
    g.tag      = string(getf(M,'tag',files(i).name));
    g.n_upd    = getf(M,'n_updates',NaN);
    g.t_med_ms = getf(M,'t_median_ms',NaN);
    g.t_p95_ms = getf(M,'t_p95_ms',NaN);
    g.t_max_ms = getf(M,'t_max_ms',NaN);
    g.miss_pct = 100*getf(M,'miss_rate',NaN);
    g.lamO2min = getf(M,'lambdaO2_min',NaN);
    g.lamO2vio = getf(M,'lambdaO2_viol',NaN);
    g.Ist_vio  = getf(M,'Ist_viol',NaN);
    g.Ib_vio   = getf(M,'Ib_viol',NaN);
    g.Vcm_vio  = getf(M,'Vcm_viol',NaN);
    g.limits_source = "current scenario (physical units)";
    g.batt_As_end = getf(M,'batt_used_As',NaN);
    g.batt_SOC_check_As=getf(M,'batt_SOC_discrepancy_As',NaN);
    g.exitflag_min=getf(M,'exitflag_min',NaN);
    g.failed_updates=getf(M,'solver_failed_updates',NaN);
    g.batt_As_max = getf(M,'batt_used_As_max',NaN);
    g.Qmax_As     = getf(M,'Qmax_As',NaN);
    g.Qmax_over   = getf(M,'Qmax_overshoot_As',NaN);
    g.H2_g        = getf(M,'H2_g',getf(M,'H2_raw_end',NaN));
    g.trk_RMSE_W  = getf(M,'track_RMSE_W',NaN);
    if isempty(rows), rows = g; else, rows(end+1) = g; end %#ok<SAGROW>
end

T = struct2table(rows);
T = sortrows(T, {'method','tag'});
disp(T);
csv = fullfile(resroot,'summary.csv');
writetable(T, csv);
fprintf('\nwrote %s\n', csv);

function v = getf(s, f, d)
if isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
