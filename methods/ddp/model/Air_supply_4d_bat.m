function dx = Air_supply_4d_bat(x,u)
R_= 8.314;  %Universal gas constant
patm = 1.01325e5;  %Atomospheric pressure
psat = 4.736e4;  %Saturation pressure 353K 
Tatm = 298.15;   %Atmosphric temperature
Tst = 353.15;    %Temperature of fuel-cell stack 
Mo2 = 32e-3; %Molar mass of oxygen
Mn2 = 28e-3; %Molar mass of nitrogen
Mv = 18.02e-3; %Molar mass of vapour
Maatm = 28.97e-3; %Molar mass of atmospheric air
gam = 1.4;  %Air specific heat ratio
Cp = 1004;   %Air specific ratio
F = 96487;  %Faraday's constant
kt = 0.0153;    %motor constant
kv = 0.0153;    %motor constant
Rcm = 0.82;    %motor constant
Jcp = 5e-5;     %compressor and motor interia
eata_cp = 0.80; %Compressor efficiency
eata_cm = 0.98; %Compressor motor mechanical efficiency
n = 381;    %Number of cells in fuel-cell stack
Afc = 280;   %Fuel-cell active area (cm^2)
Vsm = 0.02;     %Supply manifold volume
Vca = 0.01;     %Single stack cathode volume
kcain = 0.3629e-5; %Cathode inlet orifice
tm = 1.275e-4;
Cd = 0.0124;    %Cathode output throttle discharge cofficient
At = 0.002;     %Cathode output throttle area
phi_atm = 0.5;  %Average ambient air relative humidity
yo2atm = 0.21;  %Oxygen molar ratio at cathode inlet
xo2atm =yo2atm*Mo2/(yo2atm*Mo2+(1-yo2atm)*Mn2); %Oxygen mass fraction at cathode inlet
watm =Mv/(yo2atm*Mo2+(1-yo2atm)*Mn2)*(phi_atm*psat/(patm-phi_atm*psat)); %Humidity ratio of inlet air
Kapper = sqrt(mean([Mo2^2,Mn2^2,Mv^2]));
Vcm_u = 300;

c1 = (R_*Tst*kcain)/(Mo2*Vca)*(xo2atm/(1+watm));
c2 = psat;
c3 = R_*Tst/Vca;
c4 = Mo2;
c5 = Mn2;
c6 = Mv*psat;
c7 = (R_*Tst*n)/(4*Vca*F);
c8 = (R_*Tst*kcain)/(Mn2*Vca)*((1-xo2atm)/(1+watm));
c9 = eata_cm*kt*kv/(Jcp*Rcm);
c10 = Cp*Tatm/(Jcp*eata_cp);
c11 = patm;
c12 = (gam-1/gam);
c13 = eata_cm*kt/(Jcp*Rcm);
c14 = R_*Tatm/(Maatm*Vsm);
c15 = 1/eata_cp;
c16 = kcain;
c17 = (Cd*At)/sqrt(R_*Tst)*sqrt(2*gam/(gam-1));
c18 = 1/gam;
c19 = (2/(gam+1))^(gam/(gam-1));
c20 = (Cd*At)/sqrt(R_*Tst)*gam^(1/2)*(2/(gam+1))^((gam+1)/2*(gam-1));
c21 = 1/Rcm;
c22 = kv;
c23 = kcain*xo2atm/(1+watm);
c24 = n*Mo2/(4*F);

Cb = 3060;
Rsd = 10e9;
Rs = @(SOC) 0.3208*exp(-29.14*SOC) +0.04669;
Cs = @(SOC) -752.9*exp(-13.51*SOC) + 703.6;
Rf = @(SOC) 6.603*exp(-155.2*SOC) +0.04984;
Cf = @(SOC) -6056*exp(-27.12*SOC) + 4475;

x = x.*[patm;patm;1e5;patm;1;1;1];
u = u.*[300;2.2*Afc;36];

if c11/(x(1)+x(2)+c2) > c19
    % if x(4)<0
    %     x(4) = 0;
    % end
dx = [
    c1*(-x(1)-x(2)+x(4)-c2)-(c3*x(1)*psi_case1(x,c2,c11,c12,c17,c18,c19,c20))/(c4*x(1)+c5*x(2)+c6)-c7*u(2);
    c8*(-x(1)-x(2)+x(4)-c2)-(c3*x(2)*psi_case1(x,c2,c11,c12,c17,c18,c19,c20))/(c4*x(1)+c5*x(2)+c6);
    -c9*x(3)-(c10./x(3))*((x(4)./c11).^c12-1)*Weight_flow(x(3),x(4))+c13*u(1);
    c14*(1+c15*((x(4)./c11)^c12-1))*(Weight_flow(x(3),x(4))-c16*(-x(1)-x(2)+x(4)-c2));
    -x(5)/(Rsd*Cb)-u(3)/Cb/15;
    -x(6)/(Rs(x(1)*100)*Cs(x(1)*100))+u(3)/Cs(x(1)*100)/15;
    -x(7)/(Rf(x(1)*100)*Cf(x(1)*100))+u(3)/Cf(x(1)*100)/15];
else
    % if x(4)<0
    %     x(4) = 0;
    % end
dx =  [
    c1*(-x(1)-x(2)+x(4)-c2)-(c3*x(1)*psi_case2(x,c2,c11,c12,c17,c18,c19,c20))/(c4*x(1)+c5*x(2)+c6)-c7*u(2);
    c8*(-x(1)-x(2)+x(4)-c2)-(c3*x(2)*psi_case2(x,c2,c11,c12,c17,c18,c19,c20))/(c4*x(1)+c5*x(2)+c6);
    -c9*x(3)-(c10./x(3))*((x(4)./c11).^c12-1)*Weight_flow(x(3),x(4))+c13*u(1);
    c14*(1+c15*((x(4)./c11)^c12-1))*(Weight_flow(x(3),x(4))-c16*(-x(1)-x(2)+x(4)-c2))
    -x(5)/(Rsd*Cb)-u(3)/Cb/15;
    -x(6)/(Rs(x(1)*100)*Cs(x(1)*100))+u(3)/Cs(x(1)*100)/15;
    -x(7)/(Rf(x(1)*100)*Cf(x(1)*100))+u(3)/Cf(x(1)*100)/15];
end
dx = dx./[patm;patm;1e5;patm;1;1;1];



end

function psi = psi_case1(x,c2,c11,c12,c17,c18,c19,c20)
% if c11-(x(1)+x(2)+c2) > 0
%  x(1) = c11-(x(2)+c2)+1e-10;
% end
psi = c17*(x(1)+x(2)+c2)*(c11/(x(1)+x(2)+c2))^c18*sqrt(1-(c11/(x(1)+x(2)+c2))^c12);
end

function psi = psi_case2(x,c2,c11,c12,c17,c18,c19,c20)
psi = c20*(x(1)+x(2)+c2);
end
