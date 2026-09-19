function [time,Pd_base] = build_demand(scn,seg,equilibrium_power)
%BUILD_DEMAND Construct a demand trace for any benchmark scenario.
if isfield(seg,'sim_t1') && ~isempty(seg.sim_t1)
    tend=seg.sim_t1-seg.t0;
else
    tend=seg.t1-seg.t0;
end
time=(0:scn.dt:tend).';
kind="udds";
if isfield(scn,'demand') && isfield(scn.demand,'type'), kind=string(scn.demand.type); end
switch lower(kind)
    case "udds"
        U=readtable(scn.udds_file,'VariableNamingRule','preserve');
        traction=interp1(U{:,1},U{:,end},seg.t0+time,'linear','extrap');
        traction(isnan(traction))=0;
        Pd_base=equilibrium_power+scn.traction_scale*traction;
    case "step"
        Pd_base=repmat(equilibrium_power,numel(time),1);
        for j=1:numel(scn.demand.step_times_s)
            active=time>=scn.demand.step_times_s(j)-1e-12;
            Pd_base(active)=Pd_base(active)+scn.demand.step_increments_W(j);
        end
    otherwise
        error('Unknown demand type: %s',kind);
end
end
