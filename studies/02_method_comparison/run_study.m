function files=run_study()
%RUN_STUDY Stage 02: run the matched-battery-use method comparison.
here=fileparts(mfilename('fullpath')); expdir=fileparts(fileparts(here));
addpath(here,fullfile(expdir,'core'));
S=step_demand_scenario(8000,'02_method_comparison');
S.lowpath=struct('Qs',0.15,'R',100,'regenerate_functions',1);
files=cell(4,1);
files{1}=run_lowpath(S,1,16);
Mlp=eval_run(files{1});
S.segments.Qmax_As=Mlp.batt_used_As;
S.cost=struct('Qs',100,'R',diag([1 1 .01]),'regenerate_functions',1);
S.segments.name='step_demand_8kW_Qmax60p77';
files{2}=run_fmincon(S,1,10);
S.cost.regenerate_functions=0;
S.segments.name='step_demand_8kW_Qs100_Qmax60p77';
files{3}=run_ddp(S,1,50,'unif10');
files{4}=run_ddp(S,1,20,'nonunif5');
summarize_method_comparison;
end
