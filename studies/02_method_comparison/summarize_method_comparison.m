function T=summarize_method_comparison()
%SUMMARIZE_METHOD_COMPARISON Paper-facing status of the matched-use 8 kW runs.
here=fileparts(mfilename('fullpath')); expdir=fileparts(fileparts(here));
out=fullfile(expdir,'results','02_method_comparison');
addpath(fullfile(expdir,'core'));
files={fullfile(out,'lowpath','step_demand_wn16_migrated.mat'); ...
       fullfile(out,'fmincon','step_demand_8kW_Qmax60p77_unifNp10.mat'); ...
       fullfile(out,'ddp','step_demand_8kW_Qs100_Qmax60p77_unif10_N50.mat'); ...
       fullfile(out,'ddp','step_demand_8kW_Qs100_Qmax60p77_nonunif5_N20.mat')};
method={'Low-pass, omega_n=16';'fmincon uniform N=10'; ...
        'DDP uniform N=10';'DDP nonuniform N=5 (20 iter)'};
status={'benchmark; lambda constraint not enforced';'reference; one failed update'; ...
        'near-feasible';'exploratory; optimization not converged'};
n=numel(files); A=nan(n,12); sig=cell(n,1);
for i=1:n
    M=eval_run(files{i}); R=load(files{i}); so=R.simout;
    p=get_delivered(so);
    sig{i}=struct('p',p,'demand',R.demand);
    failed=0; if isfield(M,'solver_failed_updates'), failed=M.solver_failed_updates; end
    A(i,:)=[M.H2_g M.batt_used_As M.Qmax_overshoot_As M.track_RMSE_W M.lambdaO2_min ...
        M.lambdaO2_viol field_or_nan(M,'Ib_abs_max') field_or_nan(M,'Ib_viol') ...
        field_or_nan(M,'t_median_ms') field_or_nan(M,'t_p95_ms') field_or_nan(M,'miss_rate') failed];
end
T=table(method,status,A(:,1),A(:,2),A(:,3),A(:,4),A(:,5),A(:,6),A(:,7),A(:,8),A(:,9),A(:,10),A(:,11),A(:,12), ...
 'VariableNames',{'Method','Status','H2_g','Battery_As','QmaxOvershoot_As','TrackingRMSE_W','LambdaO2Min','LambdaO2Violation','IbMax_A','IbViolation_A','Median_ms','P95_ms','DeadlineMissRate','FailedUpdates'});
writetable(T,fullfile(out,'method_comparison_8kw.csv')); disp(T);
C=lines(n); f=figure('Visible','off','Color','w','Position',[100 100 900 620]); tl=tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
ax=nexttile(tl);hold(ax,'on');grid(ax,'on');plot(ax,sig{1}.demand.Time,sig{1}.demand.Data,'k--','LineWidth',1.4);
for i=1:n,plot(ax,sig{i}.p.Time,sig{i}.p.Data,'Color',C(i,:),'LineWidth',1.2);end
xlim(ax,[0 5]);ylabel(ax,'Delivered power [W]');legend(ax,[{'Demand'};method],'Location','southeast');
ax=nexttile(tl);hold(ax,'on');grid(ax,'on');
for i=1:n,d=interp1(sig{i}.demand.Time,sig{i}.demand.Data,sig{i}.p.Time,'previous','extrap');plot(ax,sig{i}.p.Time,sig{i}.p.Data(:)-d(:),'Color',C(i,:));end
xlim(ax,[0 5]);xlabel(ax,'Time [s]');ylabel(ax,'Tracking error [W]');yline(ax,0,'k:');
exportgraphics(f,fullfile(out,'method_tracking_8kw.png'),'Resolution',200);close(f);
f=figure('Visible','off','Color','w','Position',[100 100 800 430]);b=bar([A(:,9),A(:,10)]);grid on;hold on;yline(50,'r--','50 ms deadline','LineWidth',1.3);
xticklabels(method);xtickangle(15);ylabel('MPC computation time [ms]');legend(b,{'Median','95th percentile'},'Location','northwest');
exportgraphics(f,fullfile(out,'method_timing_8kw.png'),'Resolution',200);close(f);
end

function v=field_or_nan(S,name)
if isfield(S,name), v=S.(name); else, v=NaN; end
end

function p=get_delivered(so)
if isstruct(so)
    if isfield(so,'Pd'), p=so.Pd; return; end
    fc=so.Pfc; cm=so.Pcm; b=so.Pb;
else
    try, p=so.get('Pd'); if ~isempty(p), return; end, catch, end
    fc=so.get('Pfc'); cm=so.get('Pcm'); b=so.get('Pb');
end
t=fc.Time(:);
y=fc.Data(:)-interp1(cm.Time,cm.Data(:),t,'linear','extrap') ...
    +interp1(b.Time,b.Data(:),t,'linear','extrap');
p=timeseries(y,t);
end
