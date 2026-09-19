function files = run_study()
%RUN_STUDY Run the Sec. V-A Qmax sweep with the shared DDP implementation.
here=fileparts(mfilename('fullpath')); expdir=fileparts(fileparts(here));
addpath(here,fullfile(expdir,'studies','udds'),fullfile(expdir,'core'));
S=step_demand_scenario(); outdir=fullfile(expdir,'results','step_demand');
if ~isfolder(outdir), mkdir(outdir); end
files=cell(numel(S.Qmax_As),1);
for i=1:numel(S.Qmax_As)
    S.segments.Qmax_As=S.Qmax_As(i);
    S.segments.name=sprintf('step_demand_Qmax_%gAs',S.Qmax_As(i));
    files{i}=run_ddp(S,1,100,'unif10');
end
save(fullfile(outdir,'run_index.mat'),'files','S');
plot_results(outdir);
end
