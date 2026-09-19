function T = check_solver_snapshot(case_name)
% Identical state, forecast, nominal input and budget for both solvers.
if nargin<1, case_name='free'; end
assignin('base','CHECK_case',case_name);
expdir=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(expdir,'studies','03_udds'),fullfile(expdir,'core'));
assignin('base','CHECK_expdir',expdir); assignin('base','CHECK_here',pwd);
assignin('base','CHECK_path',path); assignin('base','BM_scn',udds_scenario());
mdir=fullfile(expdir,'methods','ddp');
cd(mdir); addpath(mdir,fullfile(mdir,'model'),fullfile(mdir,'functions'));
setup_ddp;
expdir=evalin('base','CHECK_expdir');
cleanup=onCleanup(@() restore()); %#ok<NASGU>
addpath(fullfile(expdir,'methods','fmincon'));
params.MPC.N=10; params.MPC.dt=.05; params.tsp_str=.05*ones(1,10);
params.J_tsp=ones(1,10); params.use_rk4=1; params.integ_h=.025;
params.ls_backtrack=1;
x0=[x_normalized;0]; Uref=repmat(u_normalized,1,10);
pd=params.Pd+linspace(0,3000,11).'; Smax=12/Ref.Ib;
case_name=evalin('base','CHECK_case');
if strcmp(case_name,'budget_active'), Smax=.1/Ref.Ib; end
out=fullfile(expdir,'results','solver_validation',case_name); if ~isfolder(out), mkdir(out); end

% Nominal and variational integrators must describe the same discrete map.
probe=Uref; probe(3,:)=.1;
[~,~,X]=evaluate_ocp(probe,x0,Uref,pd,params,Ref,Smax);
Xddp=forwardpass_exp(x0,Uref,probe-Uref,params.tsp_str,10,params,Ref);
map_error=max(abs(X-Xddp),[],'all');
assert(map_error<1e-8,'DDP/fmincon discrete maps differ by %g',map_error);
fprintf('Shared discrete map PASS: max difference %.3g\n',map_error);

rows=struct([]); solutions=struct();
for n=[50 200 500]
    params.N_iter=n;
    timer=tic;
    [~,~,~,~,Us,~,~,~,~]=solve_snapshot(x0,Uref,pd,params,Ref,Smax,Uref,zeros(15,11));
    elapsed=toc(timer);
    [J,c,X]=evaluate_ocp(Us,x0,Uref,pd,params,Ref,Smax);
    tag=sprintf('ddp_%d',n); solutions.(tag)=Us;
    g=struct('solver',string(tag),'objective',J,'max_violation',max([0;c]), ...
        'pred_battery_As',X(8,end)*Ref.Ib,'wall_s',elapsed, ...
        'exitflag',NaN,'firstorderopt',NaN);
    if isempty(rows), rows=g; else, rows(end+1)=g; end %#ok<AGROW>
    disp(g);
    save(fullfile(out,'progress.mat'),'rows','solutions','params','Ref','x0','Uref','pd','Smax','map_error');
end
opt=optimoptions('fmincon','Algorithm','sqp','Display','off','MaxIterations',500, ...
    'MaxFunctionEvaluations',30000,'OptimalityTolerance',1e-8,'ConstraintTolerance',1e-8,'StepTolerance',1e-10);
lb=repmat([0;0;-1],10,1); ub=ones(30,1);
for start={'nominal','ddp_500'}
    if strcmp(start{1},'nominal'), initial=Uref; else, initial=solutions.ddp_500; end
    timer=tic;
    [z,J,ef,info]=fmincon(@(z) evaluate_ocp(z,x0,Uref,pd,params,Ref,Smax),initial(:), ...
        [],[],[],[],lb,ub,@(z) con(z,x0,Uref,pd,params,Ref,Smax),opt);
    elapsed=toc(timer); [~,c,X]=evaluate_ocp(z,x0,Uref,pd,params,Ref,Smax);
    tag=['fmincon_from_' start{1}]; solutions.(tag)=reshape(z,3,10);
    g=struct('solver',string(tag),'objective',J,'max_violation',max([0;c]), ...
        'pred_battery_As',X(8,end)*Ref.Ib,'wall_s',elapsed, ...
        'exitflag',ef,'firstorderopt',info.firstorderopt);
    rows(end+1)=g; disp(g); %#ok<AGROW>
end
T=struct2table(rows); writetable(T,fullfile(out,'comparison.csv'));
save(fullfile(out,'comparison.mat'),'T','solutions','params','Ref','x0','Uref','pd','Smax','map_error');
end
function [c,ceq]=con(z,x,Uref,pd,params,Ref,Smax)
[~,c]=evaluate_ocp(z,x,Uref,pd,params,Ref,Smax); ceq=[];
end
function restore()
cd(evalin('base','CHECK_here')); path(evalin('base','CHECK_path'));
end
