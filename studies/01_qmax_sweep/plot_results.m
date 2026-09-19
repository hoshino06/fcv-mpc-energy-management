function T = plot_results(outdir)
%PLOT_RESULTS Plot and summarize the shared-method step-demand sweep.
Ridx=load(fullfile(outdir,'run_index.mat')); S=Ridx.S;
n=numel(Ridx.files); C=lines(n);
names=arrayfun(@(q)sprintf('Q_{max} = %g As',q),S.Qmax_As,'UniformOutput',false);
sig=cell(n,1); H2=nan(n,1); qend=nan(n,1); qmax=nan(n,1); rmse=nan(n,1); runtime=nan(n,1);
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
    rmse(i)=sqrt(mean((sig{i}.Pd.Data(:)-d(:)).^2)); runtime(i)=R.meta.whole_sim_s;
end
Qmax_As=S.Qmax_As(:); H2_g=H2; qdis_end_As=qend; qdis_max_As=qmax;
power_rmse_W=rmse; whole_sim_s=runtime;
T=table(Qmax_As,H2_g,qdis_end_As,qdis_max_As,power_rmse_W,whole_sim_s);
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
f=figure('Visible','off','Color','w','Position',[100 100 900 600]);
tiledlayout(2,1,'TileSpacing','compact'); nexttile; hold on; grid on;
plot(sig{1}.demand.Time,sig{1}.demand.Data,'k--','LineWidth',1.5);
for i=1:n, plot(sig{i}.Pd.Time,sig{i}.Pd.Data,'Color',C(i,:)); end
xlim([0 S.figure_time_s]); ylabel('Power [W]'); legend([{'Demand'},names],'Location','best');
nexttile; hold on; grid on;
for i=1:n
    ref=interp1(sig{i}.demand.Time,sig{i}.demand.Data,sig{i}.Pd.Time,'previous','extrap');
    plot(sig{i}.Pd.Time,sig{i}.Pd.Data(:)-ref(:),'Color',C(i,:));
end
xlim([0 S.figure_time_s]); xlabel('Time [s]'); ylabel('Tracking error [W]');
exportgraphics(f,fullfile(outdir,'power_tracking.png'),'Resolution',200); close(f);
f=figure('Visible','off','Color','w'); hold on; grid on;
for i=1:n, plot(sig{i}.S.Time,sig{i}.S.Data*36,'Color',C(i,:)); end
xlim([0 S.figure_time_s]); xlabel('Time [s]'); ylabel('q_{dis} [As]'); legend(names,'Location','best');
exportgraphics(f,fullfile(outdir,'fig5_qdis.png'),'Resolution',200); close(f);
f=figure('Visible','off','Color','w'); bar(H2); grid on; xticklabels(names); xtickangle(20);
ylabel('Hydrogen consumption [g]'); exportgraphics(f,fullfile(outdir,'fig6_hydrogen.png'),'Resolution',200); close(f);
disp(T);
end
