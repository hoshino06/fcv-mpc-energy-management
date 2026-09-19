function [p,t] = net_delivered_power(R)
%NET_DELIVERED_POWER Physical net bus power [W], never a filter command.
% LowPath_P is demand + compressor power, NOT delivered net power.
[p,t] = read_signal(R,'Pd');
if ~isempty(p), return; end
[pf,tf] = read_signal(R,'Pfc');
[pb,tb] = read_signal(R,'Pb');
[pc,tc] = read_signal(R,'Pcm');
p = []; t = [];
if isempty(pf) || isempty(pb) || isempty(pc), return; end
lo = max([tf(1),tb(1),tc(1)]); hi = min([tf(end),tb(end),tc(end)]);
t = unique([tf;tb;tc]); t = t(t>=lo & t<=hi);
if isempty(t), return; end
p = align(tf,pf,t) + align(tb,pb,t) - align(tc,pc,t);
end

function y = align(t,x,q)
[t,idx] = unique(t,'last'); x = x(idx);
if numel(t)==1, y = repmat(x,size(q)); else, y = interp1(t,x,q,'linear'); end
end

function [d,t] = read_signal(R,name)
d=[]; t=[];
if isstruct(R.simout)
    if ~isfield(R.simout,name), return; end
    v=R.simout.(name);
else
    try, v=R.simout.get(name); catch, return; end
end
if isa(v,'timeseries')
    d=squeeze(v.Data); t=v.Time(:);
elseif isnumeric(v) && ~isempty(v)
    d=v; t=(0:numel(v)-1).'*R.params.MPC.dt;
else, return;
end
assert(isvector(d) && numel(d)==numel(t),'Expected scalar power signal: %s',name);
d=d(:);
end
