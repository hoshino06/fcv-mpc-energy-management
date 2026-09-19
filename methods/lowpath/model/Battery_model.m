function [dx,Pb_fun] = Battery_model(Ref,n_parallel,n_series)
if nargin<2, n_parallel=15; end
if nargin<3, n_series=12; end


%x = [Vsoc;Vs;Vf]

x_normalized = sym("x_normalized",[3 1]);
x = Ref.Batt.*x_normalized;

syms Ib_normalized dIb_normalized
Ib = Ref.Ib*Ib_normalized;
dIb = dIb_normalized*Ref.Ib;
Ib = Ib/n_parallel;
dIb = dIb/n_parallel;

Vsoc = x(1);
Vs = x(2);
Vf = x(3);

Cb = 3060;
SOC = Vsoc*100;
Rsd = 10e9;
Rs = 0.3208*exp(-29.14*SOC) +0.04669;
Cs = -752.9*exp(-13.51*SOC) + 703.6;
Rf = 6.603*exp(-155.2*SOC) +0.04984;
Cf = -6056*exp(-27.12*SOC) + 4475;



dx_sym = [-Vsoc/(Rsd*Cb)-(Ib+dIb)/Cb;
    -Vs/(Rs*Cs)+(Ib+dIb)/Cs;
    -Vf/(Rf*Cf)+(Ib+dIb)/Cf];
dx = matlabFunction(dx_sym,'Vars',{x_normalized,Ib_normalized,dIb_normalized});



Voc = -1.031*exp(-35*Vsoc)+3.685+0.2156*Vsoc-0.1178*Vsoc^2+0.3201*Vsoc^3;
R0 = 0.1562*exp(-24.37*Vsoc)+0.07446;
Vb = Voc-R0*(Ib+dIb)-Vs-Vf;
Ib = Ib*n_parallel;
dIb = dIb*n_parallel;
Pb_sym = (Ib+dIb)*Vb*n_series;
Pb_fun = matlabFunction(Pb_sym,'Vars',{x_normalized,Ib_normalized,dIb_normalized});
 

end

