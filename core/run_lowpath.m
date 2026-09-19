function outfile = run_lowpath(scenario, k, wn)
%RUN_LOWPATH  Frequency-decoupling benchmark on one scenario segment.
%
%   outfile = run_lowpath(scenario, k)          % default filter (params.w_n)
%   outfile = run_lowpath(scenario, k, wn)      % override natural frequency
%
% Model: methods/lowpath/main_lowpath.slx (compiles without legacy cache).
% A low-pass filter splits the demand into an FC (low-frequency) part and a
% battery (high-frequency) part; a small NMPC tracks the resulting oxygen
% excess ratio. NOTE: this model has no per-step timing yet -> eval_run reports
% timing as unavailable until `time_ILQR` is wired into its ILQR chart.
%
% Base-workspace simulation; setup preserves the caller workspace.

if nargin < 3, wn = []; end

expdir = fileparts(fileparts(mfilename('fullpath')));      % .../experiments
bmdir  = expdir;                                            % results root
mdir = fullfile(expdir,'methods','lowpath');
assert(isfolder(mdir), 'lowpath model dir not found: %s', mdir);

assignin('base','BM_here',   pwd);
assignin('base','BM_oldpath',path);
assignin('base','BM_scn',    scenario);
assignin('base','BM_k',      k);
assignin('base','BM_wn',     wn);
assignin('base','BM_mdir',   mdir);
assignin('base','BM_bmdir',  bmdir);

try
    cd(mdir);
    addpath(mdir, fullfile(mdir,'model'), fullfile(mdir,'functions'));
    configure_model;
    setup_lowpath;

    mdir  = evalin('base','BM_mdir');  bmdir = evalin('base','BM_bmdir');
    scn   = evalin('base','BM_scn');   k = evalin('base','BM_k');  wn = evalin('base','BM_wn');
    cd(mdir);
    seg = scn.segments(k);
    if ~isempty(wn),                 params.w_n = wn;
    elseif ~isempty(seg.lowpass_tau_s), params.w_n = 1/seg.lowpass_tau_s; end

    % ---- rebuild demand for this interval (uniform N=10, terminal hold) ----
    addpath(fullfile(bmdir,'core'));
    [time,Pd_base] = build_demand(scn,seg,params.Pd);

    Ncol = params.MPC.N;
    Pd_future = zeros(numel(time), Ncol);
    for ii = 1:Ncol
        Pd_future(:,ii) = Pd_base(min((1:numel(time))+(ii-1),numel(time)));
    end
    Data.demand     = timeseries(Pd_future, time);
    Demandseries    = Data.demand;
    simulation_time = seg.t1-seg.t0;

    vars = who;
    for iv = 1:numel(vars)
        if ~startsWith(vars{iv},'BM_'), assignin('base', vars{iv}, eval(vars{iv})); end
    end

    settle = 0; if isfield(scn,'eval_settle_s'), settle = scn.eval_settle_s; end
    meta = struct('method','lowpath','tag',sprintf('%s_wn%g_migrated',seg.name,params.w_n), ...
        'timestamp',datetime('now','TimeZone','local'),'matlab',version, ...
        'computer',computer,'model','main_lowpath','segment',seg,'w_n',params.w_n, ...
        'demand_type',demand_type(scn), 'score_win_s',[settle, seg.t1 - seg.t0]);
    meta.implementation = 'methods/lowpath; corrected RK argument order';
    load_system('main_lowpath');
    % getChecksum triggers a separate code-generation compile of this legacy
    % diagram; retain source file metadata without compiling it twice.
    meta.model_file = dir(fullfile(mdir,'main_lowpath.slx'));

    t0 = tic;
    simout = sim('main_lowpath','SimulationMode','normal','StopTime',num2str(simulation_time), ...
        'ReturnWorkspaceOutputs','on','SolverType','Fixed-step','Solver','ode4', ...
        'FixedStep',num2str(scn.plant_dt));
    meta.plant_dt_s=scn.plant_dt;
    meta.whole_sim_s = toc(t0);

    demand = timeseries(Pd_base, time);
    limits = scn.limits;
    outdir = fullfile(bmdir,'results',study_name(scn),'lowpath');
    if ~isfolder(outdir), mkdir(outdir); end
    outfile = fullfile(outdir, [meta.tag '.mat']);
    save(outfile,'simout','meta','params','Ref','seg','limits','demand','-v7.3');
    fprintf('[run_lowpath] %s  (whole sim %.1f s) -> %s\n', meta.tag, meta.whole_sim_s, outfile);
    try, close_system('main_lowpath',0); catch, end

    cd(evalin('base','BM_here')); path(evalin('base','BM_oldpath'));
catch ME
    try, cd(evalin('base','BM_here')); path(evalin('base','BM_oldpath')); catch, end
    rethrow(ME);
end
end

function s=demand_type(scn)
s='udds'; if isfield(scn,'demand')&&isfield(scn.demand,'type'), s=char(scn.demand.type); end
end

function s=study_name(scn)
s='unclassified'; if isfield(scn,'study'), s=char(scn.study); end
end
