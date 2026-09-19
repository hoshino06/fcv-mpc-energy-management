function [tkp1,ykp1]= runge_Kutta(fun,dt,tk,yk)
delta_t = 0.005; %0.0001;
N = round(dt/delta_t);
for ii =1:N
k1 = fun(tk,yk);
k2 = fun(tk+delta_t/2,yk+delta_t/2*k1);
k3 = fun(tk+delta_t/2,yk+delta_t/2*k2);
k4 = fun(tk+delta_t,yk+delta_t*k3);
yk = yk+(k1+2*k2+2*k3+k4)/6*delta_t;
tk = tk+delta_t;
end
ykp1 = yk;
tkp1 = tk;
end