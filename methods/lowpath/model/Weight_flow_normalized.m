function  Wcp = Weight_flow_normalized(Ncp,Psm,x3ref,patm)
%UNTITLED2 この関数の概要をここに記述
%   詳細説明をここに記述

Ncp = Ncp*30./pi*x3ref;
Pcp = Psm*patm;
Temp = 298.15;
Baro = 101325; %Pa
dc = 0.2286;  %the compressor diameter(m)

Ncr = Ncp./sqrt(Temp./288); %the corrected compressor speed (rpm)
gam = 1.4;   %the ratio of the specific heats of the gas at constant pressure Cp./Cv

Cp = 1004;  %the specific heat capasity of air (J*kg^-1*K^-1)
Ra = 2.869e+2;  %the air gas constant (J*kg^-1*K^-1)
rho_a = 1.23;   %the air density (kg./m^3)

Uc = (pi./60)*dc*Ncr; %the compressor blade tip speed    


M = Uc./sqrt(gam*Ra*Temp);  %the inlet Mach number  
a = [-3.69906e-5 2.70399e-4 -5.36235e-4 -4.63685e-5 2.21195e-3];   %curve fitting a4,a3,a2...
b = [1.76567 -1.34837 2.44419];
c = [-9.78755e-3 0.10581 -0.42937 0.80121 -0.68344 0.43331];

% phi_max = polyval(a,M);
phi_max = a(1)*M^4+a(2)*M^3+a(3)*M^2+a(4)*M^1+a(5);
% beta = polyval(b,M);
beta = b(1)*M^2+b(2)*M+b(3);
%psi_max = polyval(c,M);
psi_max = c(6)*M^5+c(5)*M^4+c(4)*M^3+c(3)*M^2+c(2)*M+c(1);
psi = Cp*Temp*((Pcp./Baro)^((gam-1)./gam-1))./(Uc^2./2);
phi = phi_max*(1-exp(beta*(psi./psi_max-1)));

Wcr = phi*rho_a*(pi./4)*dc^2*Uc; %the air mass flow
Wcp = Wcr*(Baro./101325)./sqrt(Temp./288);% the air mass flow
