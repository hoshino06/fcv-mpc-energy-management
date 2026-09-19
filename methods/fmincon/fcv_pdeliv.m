function P = fcv_pdeliv(x, u, params, Ref)
%FCV_PDELIV  Net electrical power delivered to the bus (Pfc - Pcm + Pb) [W].
% Matches Pd_fun_case* used in the DDP cost. du = 0 (physical input in u).
x = x(:); u = u(:); z3 = zeros(3,1);
x(1:4)=max(x(1:4),1e-12);
if params.c11/(x(1)+x(2)+params.c2/Ref.pres) > params.c19*Ref.pres
    P = Pd_fun_case1(x, u, z3);
else
    P = Pd_fun_case2(x, u, z3);
end
P = P(1);
end
