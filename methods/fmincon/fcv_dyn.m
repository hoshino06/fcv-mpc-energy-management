function dx = fcv_dyn(x, u, params, Ref)
%FCV_DYN  Continuous-time 8-state FC+battery dynamics, physical (normalised) input.
%   x : [air(4); Vsoc; Vs; Vf; qdis]
%   u : [Vcm; Ist; Ib] normalised by Ref.Vcm / Ref.Ifc / Ref.Ib
% Same case selection as backwardpass_exp.m/func. u and du always enter the
% generated dynamics as (u+du), so passing du=0 and u=u_phys is exact.
x  = x(:);  u = u(:);
x(1:4)=max(x(1:4),1e-12);
z3 = zeros(3,1);
if (params.c11/Ref.pres - (x(1)+x(2)+params.c2/Ref.pres)) > 0 || x(4) < 1
    dx = zeros(8,1);
elseif params.c11/(x(1)+x(2)+params.c2/Ref.pres) > params.c19*Ref.pres
    dx = Air_supply_4d_bat_case1(x, u, z3);
else
    dx = Air_supply_4d_bat_case2(x, u, z3);
end
dx = dx(:);
end
