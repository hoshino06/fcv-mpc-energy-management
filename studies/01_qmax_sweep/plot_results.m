function T = plot_results(outdir)
%PLOT_RESULTS Plot and summarize the shared-method step-demand sweep.
Ridx=load(fullfile(outdir,'run_index.mat')); S=Ridx.S;
n=numel(Ridx.files); C=lines(n);
names=arrayfun(@(q)sprintf('Q_{max} = %g As',q),S.Qmax_As,'UniformOutput',false);
sig=cell(n,1); H2=nan(n,1); qend=nan(n,1); qmax=nan(n,1); rmse=nan(n,1); runtime=nan(n,1);
maxerr=nan(n,1); meanerr=nan(n,1); finalerr=nan(n,1); poststep_rmse=nan(n,1);
for i=1:n
    runfile=Ridx.files{i};
    if ~isfile(runfile)
        [~,name,ext]=fileparts(runfile);
        runfile=fullfile(outdir,'ddp',[name ext]);
    end
    R=load(runfile); so=R.simout; get=@(s)so.get(s);
    sig{i}=struct('Pd',get('Pd'),'Vcm',get('Vcm'),'Ifc',get('Ifc'), ...
        'Ib',get('Ib'),'S',get('S'),'H2',get('H2'),'demand',R.demand);
    H2(i)=R.Ref.Ifc*2*sig{i}.H2.Data(end)*R.params.n/(96485.3321*2);
    q=sig{i}.S.Data(:)*R.Ref.Ib; qend(i)=q(end); qmax(i)=max(q);
    d=interp1(R.demand.Time,R.demand.Data,sig{i}.Pd.Time,'previous','extrap');
    e=sig{i}.Pd.Data(:)-d(:);
    rmse(i)=sqrt(mean(e.^2)); maxerr(i)=max(abs(e)); meanerr(i)=mean(e);
    finalerr(i)=e(end); poststep_rmse(i)=sqrt(mean(e(sig{i}.Pd.Time>=2).^2));
    runtime(i)=R.meta.whole_sim_s;
end
Qmax_As=S.Qmax_As(:); H2_g=H2; qdis_end_As=qend; qdis_max_As=qmax;
power_rmse_W=rmse; power_max_abs_error_W=maxerr; power_mean_error_W=meanerr;
power_final_error_W=finalerr; power_poststep_rmse_W=poststep_rmse;
whole_sim_s=runtime;
T=table(Qmax_As,H2_g,qdis_end_As,qdis_max_As,power_rmse_W, ...
    power_max_abs_error_W,power_mean_error_W,power_final_error_W, ...
    power_poststep_rmse_W,whole_sim_s);
writetable(T,fullfile(outdir,'summary.csv'));
f=figure('Visible','off','Color','w','Position',[100 100 900 720]);
tiledlayout(3,1,'TileSpacing','compact'); labels={'V_{cm} [V]','I_{st} [A]','I_{bat} [A]'}; fields={'Vcm','Ifc','Ib'};
for k=1:3
    nexttile; hold on; grid on;
    for i=1:n, plot(sig{i}.(fields{k}).Time,sig{i}.(fields{k}).Data,'Color',C(i,:)); end
    xlim([0 S.figure_time_s]); ylabel(labels{k});
end
xlabel('Time [s]'); legend(names,'Location','best');
exportgraphics(f,fullfile(outdir,'fig4_inputs.png'),'Resolution',200); close(f);
% Paper-facing power plot: reproduce the MATLAB layout used for Fig. 3 in
% the manuscript (four supplied-power traces followed by the red demand).
paper_colors=[0 0.4470 0.7410; 0.8500 0.3250 0.0980; ...
    0.4940 0.1840 0.5560; 0.6350 0.0780 0.1840];
f=figure('Visible','off','Color','w','Position',[100 100 1400 900]);
ax=axes(f); hold(ax,'on'); box(ax,'on');
for i=1:n
    plot(ax,sig{i}.Pd.Time,sig{i}.Pd.Data,'Color',paper_colors(i,:),'LineWidth',2);
end
plot(ax,sig{1}.demand.Time,sig{1}.demand.Data,'r','LineWidth',2);
xlim(ax,[0 S.figure_time_s]); ylim(ax,[2e4 5e4]);
xlabel(ax,'Time [s]'); ylabel(ax,'Power [W]');
legend(ax,[names(:);{'Power Demand'}],'Location','best');
set(ax,'FontSize',15);
exportgraphics(f,fullfile(outdir,'power_tracking.png'),'Resolution',200); close(f);

% Diagnostic error plot is kept separate so the paper-facing figure and the
% quantitative tracking assessment do not obscure one another.
f=figure('Visible','off','Color','w','Position',[100 100 900 420]);
ax=axes(f); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
for i=1:n
    ref=interp1(sig{i}.demand.Time,sig{i}.demand.Data,sig{i}.Pd.Time,'previous','extrap');
    plot(ax,sig{i}.Pd.Time,sig{i}.Pd.Data(:)-ref(:),'Color',paper_colors(i,:),'LineWidth',1.5);
end
xlim(ax,[0 S.figure_time_s]); xlabel(ax,'Time [s]'); ylabel(ax,'P_{sys}-P_{ref} [W]');
yline(ax,0,'k:'); legend(ax,names,'Location','best');
exportgraphics(f,fullfile(outdir,'tracking_error.png'),'Resolution',200); close(f);
f=figure('Visible','off','Color','w'); hold on; grid on;
for i=1:n, plot(sig{i}.S.Time,sig{i}.S.Data*36,'Color',C(i,:)); end
xlim([0 S.figure_time_s]); xlabel('Time [s]'); ylabel('q_{dis} [As]'); legend(names,'Location','best');
exportgraphics(f,fullfile(outdir,'fig5_qdis.png'),'Resolution',200); close(f);
f=figure('Visible','off','Color','w'); bar(H2); grid on; xticklabels(names); xtickangle(20);
ylabel('Hydrogen consumption [g]'); exportgraphics(f,fullfile(outdir,'fig6_hydrogen.png'),'Resolution',200); close(f);
disp(T);
end
