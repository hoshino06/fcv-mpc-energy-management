function run_study(phase)
%RUN_STUDY Run the UDDS comparison. Phases: lowpath, nmpc, report, all.
if nargin<1, phase='all'; end
expdir=fileparts(mfilename('fullpath'));
expdir=fileparts(fileparts(expdir));
addpath(fullfile(expdir,'core'),fullfile(expdir,'studies','udds'));
S=udds_scenario(); out=fullfile(expdir,'results','udds');
if any(strcmp(phase,{'all','lowpath'}))
    records=struct([]); selected='';
    for wn=S.comparison.filter_wn
        try
            f=run_lowpath(S,1,wn); M=eval_run(f);
            feasible=all(isfinite([M.H2_g,M.batt_used_As,M.Ib_viol,M.Ist_viol,M.Vcm_viol,M.lambdaO2_viol])) ...
                && max([M.Ib_viol,M.Ist_viol,M.Vcm_viol,M.lambdaO2_viol])<1e-6 && M.batt_used_As>0;
            records(end+1).file=f; records(end).feasible=feasible; records(end).metrics=M; %#ok<AGROW>
            if feasible && isempty(selected), selected=f; end
            fprintf('LOWPATH wn=%g feasible=%d battery=%.6g As H2=%.6g g\n',wn,feasible,M.batt_used_As,M.H2_g);
        catch ME
            fprintf('LOWPATH wn=%g FAILED: %s\n',wn,ME.message);
        end
    end
    if isempty(selected)
        S.segments(1).Qmax_As=12;
        fprintf('No feasible lowpath in declared settings: using predeclared 12 As for solver comparison only.\n');
    else
        M=eval_run(selected); S.segments(1).Qmax_As=M.batt_used_As;
    end
    save(fullfile(out,'study_config.mat'),'S','records','selected');
end
if any(strcmp(phase,{'all','nmpc'}))
    R=load(fullfile(out,'study_config.mat')); S=R.S;
    for hz=S.comparison.horizons
        for ni=S.comparison.N_iter
            try
                f=run_ddp(S,1,ni,hz{1}); disp(eval_run(f));
            catch ME
                log_failure(out,sprintf('ddp_%s_%d',hz{1},ni),ME);
            end
        end
    end
    try
        f=run_fmincon(S,1,10); disp(eval_run(f));
    catch ME
        log_failure(out,'fmincon',ME);
    end
end
if any(strcmp(phase,{'all','report'}))
    collect_table;
    plot_interval(S.segments(1).name);
    close all;
end
end
function log_failure(out,name,ME)
fprintf('%s FAILED: %s\n',name,ME.message);
fid=fopen(fullfile(out,[name '_failure.txt']),'w','n','UTF-8');
fprintf(fid,'%s',getReport(ME,'extended','hyperlinks','off')); fclose(fid);
end
