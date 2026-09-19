function files = run_study()
%RUN_STUDY Stage 01: demonstrate trajectory shaping through the Qmax sweep.
here=fileparts(mfilename('fullpath')); expdir=fileparts(fileparts(here));
addpath(here,fullfile(expdir,'core'));
S=step_demand_scenario(10000,'01_qmax_sweep');
outdir=fullfile(expdir,'results','01_qmax_sweep');
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
