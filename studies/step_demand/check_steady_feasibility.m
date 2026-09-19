function T = check_steady_feasibility(outdir)
%CHECK_STEADY_FEASIBILITY Scan the legacy equilibrium map with all constraints.
here=fileparts(mfilename('fullpath'));
expdir=fullfile(here,'..','..'); mdir=fullfile(expdir,'methods','ddp');
if nargin<1, outdir=fullfile(expdir,'results','step_demand'); end
addpath(mdir,fullfile(mdir,'model'),fullfile(mdir,'functions'));
BM_scn=struct('plant',struct('start_Vcm_idx',100,'start_Ifc_idx',100,'start_SOC',.5)); %#ok<NASGU>
assignin('base','BM_scn',BM_scn); assignin('base','BM_diag_old',pwd); assignin('base','BM_diag_outdir',outdir);
cd(mdir); setup_ddp;
outdir=evalin('base','BM_diag_outdir');
target=23419.767123+[0 10000 20000];
bestErr=inf(size(target)); bestP=nan(size(target)); bestV=nan(size(target)); bestI=nan(size(target)); bestZ=nan(size(target));
maxP=-inf; maxV=nan; maxI=nan; maxZ=nan; nfeas=0;
for vc=1:300
    for ic=1:616
        xp=equil(4*(vc-1)+(1:4),ic);
        if any(~isfinite(xp)) || any(xp<=0), continue; end
        xn=[xp; .5;0;0;params.Smax]; un=[vc/Ref.Vcm;ic/Ref.Ifc;0];
        if (xn(1)/.1173+params.psat/Ref.pres)<2
            p=Pd_fun_case1(xn,un,[0;0;0]); c=constraints_case1(xn,un,[0;0;0],p,params.Smax); z=z_fun(xn(1:4),un(2),0);
        else
            p=Pd_fun_case2(xn,un,[0;0;0]); c=constraints_case2(xn,un,[0;0;0],p,params.Smax); z=z_fun(xn(1:4),un(2),0);
        end
        if ~isreal(p)||~isfinite(p)||~isreal(z)||~isfinite(z)||any(~isfinite(c))||any(c>1e-8), continue; end
        nfeas=nfeas+1;
        if p>maxP, maxP=p; maxV=vc; maxI=ic; maxZ=z; end
        for j=1:numel(target)
            e=abs(target(j)-p);
            if e<bestErr(j), bestErr(j)=e; bestP(j)=p; bestV(j)=vc; bestI(j)=ic; bestZ(j)=z; end
        end
    end
end
Demand_W=target(:); ClosestFeasible_W=bestP(:); Error_W=bestP(:)-target(:);
VcmIndex=bestV(:); IfcIndex=bestI(:); LambdaO2=bestZ(:);
T=table(Demand_W,ClosestFeasible_W,Error_W,VcmIndex,IfcIndex,LambdaO2);
writetable(T,fullfile(outdir,'steady_feasibility.csv'));
fprintf('Feasible grid points: %d; max %.3f W at (%g,%g), lambda=%.4f\n',nfeas,maxP,maxV,maxI,maxZ); disp(T);
cd(evalin('base','BM_diag_old'));
end
