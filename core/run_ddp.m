function outfile = run_ddp(scenario, k, N_iter, horizon)
%RUN_DDP  DDP-based NMPC on one scenario interval. One driver for every horizon.
%
%   outfile = run_ddp(scenario, k, N_iter [, horizon])
%     scenario : struct from a studies/<name> scenario function
%     k        : index into scenario.segments
%     N_iter   : augmented-Lagrangian iteration cap
%     horizon  : 'nonunif5' (default) | 'unif10' | 'unif5' | a 1xN step vector [s]
%
%   Control (sampling) period is ALWAYS scenario.dt (= 0.05 s) for every
%   horizon; `horizon` only sets the PREDICTION step schedule tsp_str:
%     'nonunif5' : exponential N=5, first step 0.05 s, 0.5 s span  (proposed)
%     'unif10'   : uniform  N=10 x 0.05 s, 0.5 s span              (baseline)
%     'unif5'    : uniform  N=5  x 0.05 s, 0.25 s span             (short-horizon)
%
%   Model: methods/ddp/main_ddp_N<N>.slx (chart hard-coded to that N).
%   NOTE: setup_ddp.m preserves BM_scn across its workspace reset. Other runner
%   inputs remain stashed in base. Solver knobs come from scenario.solver.

if nargin < 4 || isempty(horizon), horizon = 'nonunif5'; end

benchdir = fileparts(mfilename('fullpath')); % .../fcv-mpc-energy-management/core
expdir   = fileparts(benchdir);              % .../fcv-mpc-energy-management
mdir     = fullfile(expdir,'methods','ddp');
assert(isfolder(mdir), 'DDP model dir not found: %s', mdir);

assignin('base','BM_here',pwd);        assignin('base','BM_oldpath',path);
assignin('base','BM_scn',scenario);    assignin('base','BM_k',k);
assignin('base','BM_Niter',N_iter);    assignin('base','BM_horizon',horizon);
assignin('base','BM_mdir',mdir);       assignin('base','BM_expdir',expdir);

try
    cd(mdir);
    addpath(mdir, fullfile(mdir,'model'), fullfile(mdir,'functions'));
    setup_ddp;                                   % <- `clear`; sets params, Ref, tsp_str1, ...

    mdir=evalin('base','BM_mdir'); expdir=evalin('base','BM_expdir');
    scn=evalin('base','BM_scn'); k=evalin('base','BM_k');
    N_iter=evalin('base','BM_Niter'); horizon=evalin('base','BM_horizon');
    cd(mdir);
    seg = scn.segments(k);

    % ---- resolve the prediction step schedule ----
    if ischar(horizon) || isstring(horizon)
        switch lower(string(horizon))
            case "nonunif5", tsp = params.tsp_str(:).';  hl = 'nonunif5';
            case "unif10",   tsp = 0.05*ones(1,10);      hl = 'unif10';
            case "unif5",    tsp = 0.05*ones(1,5);       hl = 'unif5';
            otherwise, error('unknown horizon "%s"', horizon);
        end
    else
        tsp = horizon(:).';
        if max(abs(diff(tsp))) < 1e-9, hl = sprintf('unif%d', numel(tsp));
        else, hl = sprintf('nonunif%d', numel(tsp)); end
    end
    Nh = numel(tsp);
    model = sprintf('main_ddp_N%d', Nh);
    assert(isfile(fullfile(mdir,[model '.slx'])), 'no DDP model for N=%d (%s.slx)', Nh, model);

    % ---- overrides ----
    params.use_rk4      = scn.solver.use_rk4;
    params.integ_h      = scn.solver.integ_h;
    params.ls_backtrack = scn.solver.ls_backtrack;
    params.N_iter       = N_iter;
    params.Smax         = seg.Qmax_As / Ref.Ib;
    params.MPC.N        = Nh;
    params.MPC.dt       = scn.dt;                 % control period, always 0.05 s
    params.tsp_str      = tsp;
    params.J_tsp        = tsp / 0.05;
    tsp_str1            = cumsum(tsp);

    % ---- interval demand (control period = scn.dt; lookahead uses tsp) ----
    addpath(fullfile(expdir,'core'));
    [time,Pd_base] = build_demand(scn,seg,params.Pd);

    lookahead = [0 tsp_str1];
    Pd_future = zeros(numel(time), Nh+1);
    for ii = 1:Nh+1
        Pd_future(:,ii) = interp1(time, Pd_base, min(time + lookahead(ii),time(end)), 'linear');
    end
    Data.demand = timeseries(Pd_future, time);
    Demandseries = Data.demand;
    simulation_time = seg.t1-seg.t0; % retain preview data, stop at the scored endpoint

    vars = who;
    for iv = 1:numel(vars)
        if ~startsWith(vars{iv},'BM_'), assignin('base', vars{iv}, eval(vars{iv})); end
    end

    settle = 0; if isfield(scn,'eval_settle_s'), settle = scn.eval_settle_s; end
    meta = struct('method','ddp','tag',sprintf('%s_%s_N%d',seg.name,hl,N_iter), ...
        'timestamp',datetime('now','TimeZone','local'),'matlab',version,'computer',computer, ...
        'model',model,'segment',seg,'solver',scn.solver,'N_iter',N_iter,'horizon',hl, ...
        'tsp_str',tsp,'ctrl_dt_s',scn.dt,'demand_type',demand_type(scn), ...
        'score_win_s',[settle, seg.t1 - seg.t0]);
    load_system(model);
    meta.model_file = dir(fullfile(mdir,[model '.slx']));

    t0 = tic;
    simout = sim(model,'StopTime',num2str(simulation_time),'ReturnWorkspaceOutputs','on', ...
        'SolverType','Fixed-step','Solver','ode4','FixedStep',num2str(scn.plant_dt));
    meta.plant_dt_s=scn.plant_dt;
    meta.whole_sim_s = toc(t0);

    demand = timeseries(Pd_base, time);
    limits = scn.limits;
    outdir = fullfile(expdir,'results',study_name(scn),'ddp');
    if ~isfolder(outdir), mkdir(outdir); end
    outfile = fullfile(outdir, [meta.tag '.mat']);
    save(outfile,'simout','meta','params','Ref','seg','limits','demand','-v7.3');
    fprintf('[run_ddp] %s  (%s, whole sim %.1f s) -> %s\n', meta.tag, model, meta.whole_sim_s, outfile);
    try, close_system(model,0); catch, end

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
