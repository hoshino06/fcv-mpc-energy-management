function outfile = run_fmincon(scenario, k, steps)
%RUN_FMINCON  Reference NMPC on one interval via a hand-rolled direct-shooting
%             fmincon, with an explicit (possibly non-uniform) step schedule.
%
%   outfile = run_fmincon(scenario, k [, steps])
%     steps : 1xN prediction step widths [s]  (default = uniform 10 x 0.05 s), OR
%             a scalar integer N -> uniform N x 0.05 s, OR
%             'nonuniform' -> the setting's params.tsp_str (non-uniform N=5,
%             same 0.5 s span, same 0.05 s FIRST step). This is the key contrast:
%             non-uniform N=5 vs uniform N=10 at matched first-step + horizon.
%
% Same OCP as the proposed DDP: decision var = input trajectory U (nu x N),
% Generated DDP stage/terminal costs and all 15 constraints are evaluated
% by evaluate_ocp. du is relative to the previous nominal input trajectory.
% The first optimized input is held for one 50 ms period.
%
% Saves an eval_run-compatible struct to results/fmincon/<seg>_<tag>.mat.

if nargin < 3 || isempty(steps), steps = 0.05*ones(1,10); end

expdir = fileparts(fileparts(mfilename('fullpath')));      % .../experiments
bmdir  = expdir;
mdir   = fullfile(expdir,'methods','ddp');
assignin('base','BM_here',pwd); assignin('base','BM_oldpath',path);
assignin('base','BM_scn',scenario); assignin('base','BM_k',k);
assignin('base','BM_steps',steps); assignin('base','BM_mdir',mdir); assignin('base','BM_bmdir',bmdir);

try
    cd(mdir);
    addpath(mdir, fullfile(mdir,'model'), fullfile(mdir,'functions'), ...
            fullfile(bmdir,'methods','fmincon'), fullfile(bmdir,'core'));
    setup_ddp;

    mdir=evalin('base','BM_mdir'); bmdir=evalin('base','BM_bmdir');
    scn=evalin('base','BM_scn'); k=evalin('base','BM_k'); steps=evalin('base','BM_steps');
    if ischar(steps) || isstring(steps)
        assert(strcmpi(steps,'nonuniform'),'steps string must be ''nonuniform''');
        steps = params.tsp_str(:).';
    elseif isscalar(steps)
        steps = 0.05*ones(1,round(steps));
    end
    steps = steps(:).';  Nh = numel(steps);
    cd(mdir);
    seg = scn.segments(k);
    dt  = scn.dt;
    Smax = seg.Qmax_As / Ref.Ib;
    lamO2min = scn.limits.lambda_O2_min;
    params.MPC.dt=dt; params.MPC.N=Nh;
    params.tsp_str=steps; params.J_tsp=steps/dt;
    params.integ_h=scn.solver.integ_h;

    % --- interval demand ---
    addpath(fullfile(bmdir,'core'));
    [tg,Pd_base] = build_demand(scn,seg,params.Pd);
    tcum = [0 cumsum(steps)];                        % same prediction nodes as DDP

    x0 = [x_normalized(:); 0];
    u0 = u_normalized(:);
    Imax=scn.limits.Ib_abs_max;
    lb = [0; 0; -Imax/Ref.Ib];  ub = [params.Vcm_u/Ref.Vcm; 2.2*params.Afc/Ref.Ifc; Imax/Ref.Ib];
    LB = repmat(lb,Nh,1);  UB = repmat(ub,Nh,1);
    opt = optimoptions('fmincon','Display','off','Algorithm','sqp', ...
        'MaxFunctionEvaluations',6e3,'MaxIterations',300,'StepTolerance',1e-8);

    tc=tg(tg<=seg.t1-seg.t0+1e-9);
    Kc = numel(tc);
    L = struct('t',tc);
    Z = zeros(Kc,1);
    [L.Pd,L.Pfc,L.Pb,L.Pcm,L.z,L.Ifc,L.Ib,L.Vcm,L.S,L.Vsoc,L.Vs,L.Vf,L.H2rate,L.tsolve,L.exitflag] = ...
        deal(Z,Z,Z,Z,Z,Z,Z,Z,Z,Z,Z,Z,Z,Z,Z);

    x = x0(:);  Uwarm = repmat(u0,1,Nh); Uref=Uwarm;
    L.firstorderopt=Z; L.constrviolation=Z;
    for c = 1:Kc
        pd = interp1(tg,Pd_base,min(tg(c)+tcum(:),tg(end)),'linear');
        z0 = Uwarm(:);
        cost2 = @(z) evaluate_ocp(z,x,Uref,pd,params,Ref,Smax);
        conf  = @(z) shared_constraints(z,x,Uref,pd,params,Ref,Smax);
        tt = tic;
        [zsol,~,ef,info] = fmincon(cost2, z0, [],[],[],[], LB, UB, conf, opt);
        L.tsolve(c) = toc(tt);  L.exitflag(c) = ef;
        L.firstorderopt(c)=info.firstorderopt; L.constrviolation(c)=info.constrviolation;
        Usol = reshape(zsol, 3, Nh);
        mv = Usol(:,1);
        Uref=Usol;
        Uwarm = [Usol(:,2:end), Usol(:,end)];

        L.Pd(c) = fcv_pdeliv(x, mv, params, Ref);
        Pb = Pb_fun(x(5:7), mv(3), 0); Pcm = Pcm_fun(x(1:4), mv(1), 0);
        L.Pb(c)=Pb; L.Pcm(c)=Pcm; L.Pfc(c)=L.Pd(c)+Pcm-Pb;
        L.z(c) = z_fun(x(1:4), mv(2), 0);
        L.Ifc(c)=mv(2)*Ref.Ifc; L.Ib(c)=mv(3)*Ref.Ib; L.Vcm(c)=mv(1)*Ref.Vcm;
        L.S(c) = x(8); L.Vsoc(c)=x(5); L.Vs(c)=x(6); L.Vf(c)=x(7); L.H2rate(c)=mv(2);
        x = rk4_step(@(xx) fcv_dyn(xx,mv,params,Ref), x, dt, max(1,ceil(dt/scn.plant_dt)));
    end
    H2acc = [0;cumsum(diff(tc).*L.H2rate(1:end-1))];

    mk = @(d) timeseries(d, tc);
    simout = struct('time_ILQR',mk(L.tsolve),'Pd',mk(L.Pd),'Pfc',mk(L.Pfc),'Pb',mk(L.Pb), ...
        'Pcm',mk(L.Pcm),'z',mk(L.z),'Ifc',mk(L.Ifc),'Ib',mk(L.Ib),'Vcm',mk(L.Vcm),'S',mk(L.S),'H2',mk(H2acc), ...
        'Battery_state',timeseries([L.Vsoc L.Vs L.Vf], tc), ...
        'exitflag',mk(L.exitflag),'firstorderopt',mk(L.firstorderopt), ...
        'constrviolation',mk(L.constrviolation));

    settle = 0; if isfield(scn,'eval_settle_s'), settle = scn.eval_settle_s; end
    if abs(max(diff(steps))) < 1e-9, hl = sprintf('unifNp%d', Nh);
    else, hl = sprintf('nonunifN%d', Nh); end
    meta = struct('method','fmincon','tag',sprintf('%s_%s',seg.name,hl), ...
        'timestamp',datetime('now','TimeZone','local'),'matlab',version,'computer',computer, ...
        'model','fmincon-shooting','segment',seg,'steps',steps,'demand_type',demand_type(scn), ...
        'score_win_s',[settle, seg.t1 - seg.t0], 'exitflag_min',min(L.exitflag));
    demand = timeseries(Pd_base, tg); limits = scn.limits;
    outdir = fullfile(bmdir,'results',study_name(scn),'fmincon');
    if ~isfolder(outdir), mkdir(outdir); end
    outfile = fullfile(outdir, [meta.tag '.mat']);
    save(outfile,'simout','meta','params','Ref','seg','limits','demand','-v7.3');
    fprintf('[run_fmincon] %s  N=%d  median solve %.3f s  min exitflag %d -> %s\n', ...
        meta.tag, Nh, median(L.tsolve), min(L.exitflag), outfile);

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

% ---------------------------------------------------------------------------
function y = rk4_step(f, x, T, nsub)
h = T/nsub; y = x(:);
for s = 1:nsub
    k1=f(y); k2=f(y+0.5*h*k1); k3=f(y+0.5*h*k2); k4=f(y+h*k3);
    y = y + (h/6)*(k1+2*k2+2*k3+k4);
end
end

function [c,ceq]=shared_constraints(z,x,Uref,pd,params,Ref,Smax)
[~,c]=evaluate_ocp(z,x,Uref,pd,params,Ref,Smax); ceq=[];
end
