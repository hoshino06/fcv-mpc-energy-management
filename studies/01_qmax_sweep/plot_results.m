function T = plot_results(outdir)
%PLOT_RESULTS Plot and summarize the shared-method step-demand sweep.
Ridx=load(fullfile(outdir,'run_index.mat')); S=Ridx.S;
n=numel(Ridx.files);
names=arrayfun(@(q)sprintf('$Q_{\\max} = %g$ As',q),S.Qmax_As,'UniformOutput',false);
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
paper_colors=[0 0.4470 0.7410; 0.8500 0.3250 0.0980; ...
    0.4940 0.1840 0.5560; 0.6350 0.0780 0.1840];

% 01: manuscript Fig. 3 -- supplied power and reference demand.
f=paper_figure([720 270]); ax=axes(f); hold(ax,'on');
for i=1:n
    plot(ax,sig{i}.Pd.Time,sig{i}.Pd.Data,'Color',paper_colors(i,:),'LineWidth',1.1);
end
plot(ax,sig{1}.demand.Time,sig{1}.demand.Data,'r','LineWidth',1.1);
xlim(ax,[0 S.figure_time_s]); ylim(ax,[2e4 5e4]);
xlabel(ax,'Time [s]'); ylabel(ax,'Power [W]'); paper_axes(ax);
legend(ax,[names(:);{'Power Demand'}],'Location','northeast', ...
    'Interpreter','latex','FontName','Times New Roman','FontSize',8);
exportgraphics(f,fullfile(outdir,'01_power_tracking.png'),'Resolution',600); close(f);

% 02: manuscript Fig. 4 -- three control inputs.
f=paper_figure([720 500]); tl=tiledlayout(f,3,1,'TileSpacing','compact','Padding','compact');
labels={'Voltage [V]','Current [A]','Current [A]'}; fields={'Vcm','Ifc','Ib'};
ylims={[90 150],[90 250],[-40 40]}; yticks_={90:20:150,100:50:250,-40:20:40};
for k=1:3
    ax=nexttile(tl); hold(ax,'on');
    for i=1:n
        plot(ax,sig{i}.(fields{k}).Time,sig{i}.(fields{k}).Data, ...
            'Color',paper_colors(i,:),'LineWidth',1.0);
    end
    xlim(ax,[0 S.figure_time_s]); ylim(ax,ylims{k}); yticks(ax,yticks_{k});
    ylabel(ax,labels{k}); paper_axes(ax);
    if k<3, ax.XTickLabel=[]; else, xlabel(ax,'Time [s]'); end
    if k==1
        legend(ax,names,'Location','southeast','Interpreter','latex', ...
            'FontName','Times New Roman','FontSize',8);
    end
end
exportgraphics(f,fullfile(outdir,'02_control_inputs.png'),'Resolution',600); close(f);

% 03: manuscript Fig. 5 -- cumulative battery discharge.
f=paper_figure([720 300]); ax=axes(f); hold(ax,'on');
for i=1:n
    plot(ax,sig{i}.S.Time,sig{i}.S.Data*36,'Color',paper_colors(i,:),'LineWidth',1.1);
end
xlim(ax,[0 S.figure_time_s]); ylim(ax,[-10 80]);
xlabel(ax,'Time [s]'); ylabel(ax,'$q_{\mathrm{dis}}$ [As]','Interpreter','latex'); paper_axes(ax);
legend(ax,names,'Location','northeast','Interpreter','latex', ...
    'FontName','Times New Roman','FontSize',8);
exportgraphics(f,fullfile(outdir,'03_battery_discharge.png'),'Resolution',600); close(f);

% 04: manuscript Fig. 6 -- hydrogen consumption and relative increase.
f=paper_figure([720 270]); ax=axes(f);
b=bar(ax,H2,'FaceColor','flat','BarWidth',0.72); b.CData=paper_colors;
ylim(ax,[3 3.5]); yticks(ax,3:0.1:3.5); ylabel(ax,'Hydrogen consumption [g]');
xticks(ax,1:n); xticklabels(ax,names); ax.TickLabelInterpreter='latex';
paper_axes(ax);
relative_pct=(H2/H2(1)-1)*100;
labels_pct=arrayfun(@(v)sprintf('%+.1f\\%%',v),relative_pct,'UniformOutput',false);
labels_pct{1}='0.0\%';
text(ax,b.XEndPoints,b.YEndPoints+0.012,labels_pct,'HorizontalAlignment','center', ...
    'VerticalAlignment','bottom','Interpreter','latex','FontName','Times New Roman','FontSize',8);
exportgraphics(f,fullfile(outdir,'04_hydrogen_consumption.png'),'Resolution',600); close(f);

% Diagnostic error plot is kept separate so the paper-facing figure and the
% quantitative tracking assessment do not obscure one another.
f=paper_figure([720 300]); ax=axes(f); hold(ax,'on');
for i=1:n
    ref=interp1(sig{i}.demand.Time,sig{i}.demand.Data,sig{i}.Pd.Time,'previous','extrap');
    plot(ax,sig{i}.Pd.Time,sig{i}.Pd.Data(:)-ref(:),'Color',paper_colors(i,:),'LineWidth',1.0);
end
xlim(ax,[0 S.figure_time_s]); xlabel(ax,'Time [s]'); ylabel(ax,'P_{sys}-P_{ref} [W]');
yline(ax,0,'k:','LineWidth',0.7); paper_axes(ax);
legend(ax,names,'Location','southwest','Interpreter','latex', ...
    'FontName','Times New Roman','FontSize',8);
exportgraphics(f,fullfile(outdir,'01a_tracking_error.png'),'Resolution',600); close(f);
disp(T);
end

function f=paper_figure(sz)
f=figure('Visible','off','Color','w','Units','pixels','Position',[100 100 sz]);
end

function paper_axes(ax)
set(ax,'FontName','Times New Roman','FontSize',10,'LineWidth',0.7, ...
    'TickDir','in','Box','on','XMinorTick','off','YMinorTick','off');
end
