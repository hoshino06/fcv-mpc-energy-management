function [dx_case1,dx_case2,Pcm,Vfc,z_fun,params] = equi()
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
c14 = R_*Tatm/(Maatm*Vsm)*kcain;
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

% c9 = (R_*Tst*kcain/(Mo2*Vca))*(xo2atm/(1+watm)) + (R_*Tst*kcain/(Mn2*Vca))*((1-xo2atm)/(1+watm)) + (R_*Tst/Vca)*Cd*At/sqrt(R_*Tst)*gam^(1/2)*((2/(gam+1))^(gam+1/(2*(gam-1))))/Kapper;
% c10 = (R_*Tst*kcain/(Mo2*Vca))*(xo2atm/(1+watm))+(R_*Tst*kcain/(Mn2*Vca))*((1-xo2atm)/(1+watm));
% c11 = psat*(R_*Tst/Vca)*Cd*At/sqrt(R_*Tst)*gam^(1/2)*((2/(gam+1))^(gam+1/(2*(gam-1))))/Kapper;
% c12 = R_*Tst*n/(4*Vca*F);
% c20 = kcain*xo2atm/(1+watm);
% c21 = n*Mo2/(4*F);

variables = who;
for i = 1:size(variables)
    params.(variables{i}) = eval(variables{i});
end
dx_case1 = @(x,Vcm,Ifc) [
    c1*(-x(1)-x(2)+x(4)-c2)-(c3*x(1)*psi_case1(x,c2,c11,c12,c17,c18,c19,c20))/(c4*x(1)*c5*x(2)+c6)-c7*Ifc;
    c8*(-x(1)-x(2)+x(4)-c2)-(c3*x(2)*psi_case1(x,c2,c11,c12,c17,c18,c19,c20))/(c4*x(1)*c5*x(2)+c6);
    -c9*x(3)-(c10./x(3))*((x(4)./c11).^c12-1)*Weight_flow(x(3),x(4))+c13*+c13*Vcm;
    c14*(1+c15*((x(4)./c11)^c12-1))*(Weight_flow(x(3),x(4))-c16*(-x(1)-x(2)+x(4)-c2))];

dx_case2 = @(x,Vcm,Ifc) [
    c1*(-x(1)-x(2)+x(4)-c2)-(c3*x(1)*psi_case2(x,c2,c11,c12,c17,c18,c19,c20))/(c4*x(1)*c5*x(2)+c6)-c7*Ifc;
    c8*(-x(1)-x(2)+x(4)-c2)-(c3*x(2)*psi_case2(x,c2,c11,c12,c17,c18,c19,c20))/(c4*x(1)*c5*x(2)+c6);
    -c9*x(3)-(c10./x(3))*((x(4)./c11).^c12-1)*Weight_flow(x(3),x(4))+c13*Vcm;
    c14*(1+c15*((x(4)./c11)^c12-1))*(Weight_flow(x(3),x(4))-c16*(-x(1)-x(2)+x(4)-c2))];
x = sym("x",[4 1]);
% syms Vcm
% dx = dx_case2(x,Vcm,120);
dx = dx_case1(x,1.143851260677745e+02,120);
equation = [dx(1)==0;dx(2)==0;dx(3)==0;dx(4)==0];
% equation = subs(equation,x(3:4),[4.982436134570045e+03;1.386310872379161e+05]);
flag1 = 0;
flag2 = 0;
flag3 = 1;
count = 1;
while flag1==0 || flag2 ==0 || flag3 == 1
Answer = vpasolve(equation,x,'Random',true);
Ans = [Answer.x1;Answer.x2;Answer.x3;Answer.x4]
flag1 = isreal(Ans);
flag2 = all(Ans(:)>0);
flag3 = isempty(Ans);
if flag1==1 && flag3 ==0
Data(:,count) = double(Ans);
end
count = count+1
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Pcm = @(x,Vcm) Vcm*(Vcm-kv*x(3))/(Rcm*eata_cm);

Vfc.case1 = @(x,Ifc) Vstack1(x,Ifc,Tst,psat,patm,Tatm,Afc,tm,n);
Vfc.case2 = @(x,Ifc) Vstack2(x,Ifc,Tst,psat,patm,Tatm,Afc,tm,n);
z_fun = @(x,Ifc) c23*(x(4)-x(3))/(c24*Ifc);
end

function Vfc = Vstack1(x,Ifc,Tfc,psat,patm,Tatm,Afc,tm,n)
po2 = x(1);
pH2 = x(3)-psat;


E = 1.229-0.85e-3*(Tfc-Tatm)+4.3085e-5*Tfc*(log(pH2/patm)+1/2*log(po2/patm));
v0 = 0.279-8.5e-4*(Tfc-Tatm)+4.308e-5*Tfc*(log((x(3)-psat)/patm)+1/2*log(0.1173*(x(3)-psat)/patm));
va = (-1.68e-5*Tfc+1.618e-2)*(1.01325*(po2/0.1173+psat)/patm)^2+(1.8e-4*Tfc-0.166)*(1.01325*(po2/0.1173+psat)/patm)+(-5.8e-4*Tfc+0.5736);
a1 = 10;
I = Ifc/Afc;
vact = v0+va*(1-exp(-a1*I));
b11 = 5.139e-3;
b12 = 3.26e-3;
b2 = 350;
mu_m = 14; % vapor 100%
sig_m = (b11*mu_m-b12)*exp(b2*(1/303-1/Tfc));
Rohm = tm/sig_m;
vohm = I*Rohm;

%if 1.01325*(po2/0.1173+psat) < 2*patm
a2 = (7.16e-4*Tfc-0.622)*1.01325*(po2/0.1173+psat)/patm+(-1.45e-3*Tfc+1.68);
%else 
% a2 = (8.66e-5*Tfc-0.068)*1.01325*(po2/0.1173+psat)/patm+(-1.6e-4*Tfc+0.54);
%end

Imax = 2.2;
a3 = 2;
vconc = I*(a2*I/Imax)^a3;

vfc = E-vact-vohm-vconc;
Vfc = vfc*n;
end

function Vfc = Vstack2(x,Ifc,Tfc,psat,patm,Tatm,Afc,tm,n)
po2 = x(1);
pH2 = x(3)-psat;


E = 1.229-0.85e-3*(Tfc-Tatm)+4.3085e-5*Tfc*(log(pH2/patm)+1/2*log(po2/patm));
v0 = 0.279-8.5e-4*(Tfc-Tatm)+4.308e-5*Tfc*(log((x(3)-psat)/patm)+1/2*log(0.1173*(x(3)-psat)/patm));
va = (-1.68e-5*Tfc+1.618e-2)*(1.01325*(po2/0.1173+psat)/patm)^2+(1.8e-4*Tfc-0.166)*(1.01325*(po2/0.1173+psat)/patm)+(-5.8e-4*Tfc+0.5736);
a1 = 10;
I = Ifc/Afc;
vact = v0+va*(1-exp(-a1*I));
b11 = 5.139e-3;
b12 = 3.26e-3;
b2 = 350;
mu_m = 14; % vapor 100%
sig_m = (b11*mu_m-b12)*exp(b2*(1/303-1/Tfc));
Rohm = tm/sig_m;
vohm = I*Rohm;

%if 1.01325*(po2/0.1173+psat) < 2*patm
%a2 = (7.16e-4*Tfc-0.622)*1.01325*(po2/0.1173+psat)/patm+(-1.45e-3*Tfc+1.68);
%else 
a2 = (8.66e-5*Tfc-0.068)*1.01325*(po2/0.1173+psat)/patm+(-1.6e-4*Tfc+0.54);
%end

Imax = 2.2;
a3 = 2;
vconc = I*(a2*I/Imax)^a3;

vfc = E-vact-vohm-vconc;
Vfc = vfc*n;
end

function psi = psi_case1(x,c2,c11,c12,c17,c18,c19,c20)
psi = c17*(x(1)+x(2)+c2)*(c11/(x(1)+x(2)+c2))^c18*sqrt(1-(c11/(x(1)+x(2)+c2))^c12);
end

function psi = psi_case2(x,c2,c11,c12,c17,c18,c19,c20)
psi = c20*(x(1)+x(2)+c2);
end
