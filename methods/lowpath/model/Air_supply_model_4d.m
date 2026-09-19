function [dx_case1,dx_case2,Pcm,Pfc,z_fun,params] = Air_supply_model_4d(Ref)
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

% load value_equil_4d.mat
% 
% E(1) = (Kapper-c5)*equil(1)+c2-c6;
% E(2) = equil(1)-E(1);
% E(3) = equil(2);
% E(4) = equil(3);


variables = who;
for i = 1:size(variables,1)
    params.(variables{i}) = eval(variables{i});
end
% Ref.x = [patm;patm;1e5;patm];
% Ref.pres = patm;
% Ref.P = 1e4;
% Ref.z = 2;
% Ref.Vcm = Vcm_u;
% Ref.Ifc = 2.2*Afc;
% 
% Ref.x = Ref.x/100;
% Ref.P = Ref.P/100;
% Ref.z = Ref.z/100;
% Ref.Vcm = Ref.Vcm/100;
% Ref.Ifc = Ref.Ifc/100;

syms Vcm_normalized Ifc_normalized dVcm_normalized dIfc_normalized
x_normalized = sym("x_normalized",[4 1]);
x = Ref.x.*x_normalized;
Vcm = Vcm_normalized*Ref.Vcm;
dVcm = dVcm_normalized*Ref.Vcm;
Ifc = Ifc_normalized*Ref.Ifc;
dIfc = dIfc_normalized*Ref.Ifc; 

dx_case1_sym = [
    c1*(-x(1)-x(2)+x(4)-c2)-(c3*x(1)*psi_case1(x,c2,c11,c12,c17,c18,c19,c20))/(c4*x(1)+c5*x(2)+c6)-c7*(Ifc+dIfc);
    c8*(-x(1)-x(2)+x(4)-c2)-(c3*x(2)*psi_case1(x,c2,c11,c12,c17,c18,c19,c20))/(c4*x(1)+c5*x(2)+c6);
    -c9*x(3)-(c10./x(3))*((x(4)./c11).^c12-1)*Weight_flow(x(3),x(4))+c13*(Vcm+dVcm);
    c14*(1+c15*((x(4)./c11)^c12-1))*(Weight_flow(x(3),x(4))-c16*(-x(1)-x(2)+x(4)-c2))]./Ref.x;



dx_case2_sym = [
    c1*(-x(1)-x(2)+x(4)-c2)-(c3*x(1)*psi_case2(x,c2,c11,c12,c17,c18,c19,c20))/(c4*x(1)+c5*x(2)+c6)-c7*(Ifc+dIfc);
    c8*(-x(1)-x(2)+x(4)-c2)-(c3*x(2)*psi_case2(x,c2,c11,c12,c17,c18,c19,c20))/(c4*x(1)+c5*x(2)+c6);
    -c9*x(3)-(c10./x(3))*((x(4)./c11).^c12-1)*Weight_flow(x(3),x(4))+c13*(Vcm+dVcm);
    c14*(1+c15*((x(4)./c11)^c12-1))*(Weight_flow(x(3),x(4))-c16*(-x(1)-x(2)+x(4)-c2))]./Ref.x;


dx_case1 = matlabFunction(dx_case1_sym,'Vars',{x_normalized,Vcm_normalized,dVcm_normalized,Ifc_normalized,dIfc_normalized});
dx_case2 = matlabFunction(dx_case2_sym,'Vars',{x_normalized,Vcm_normalized,dVcm_normalized,Ifc_normalized,dIfc_normalized});

Pcm_sym = (Vcm+dVcm)*((Vcm+dVcm)-kv*x(3))/(Rcm*eata_cm);
Pcm = matlabFunction(Pcm_sym,'Vars',{x_normalized,Vcm_normalized,dVcm_normalized});


Vfc.case1_sym =  Vstack1(x,(Ifc+dIfc),Tst,psat,patm,Tatm,Afc,tm,n);
Vfc.case2_sym = Vstack2(x,(Ifc+dIfc),Tst,psat,patm,Tatm,Afc,tm,n);
Vfc.case1 = matlabFunction(Vfc.case1_sym,'Vars',{x_normalized,Ifc_normalized,dIfc_normalized});
Vfc.case2 = matlabFunction(Vfc.case2_sym,'Vars',{x_normalized,Ifc_normalized,dIfc_normalized});
Pfc_case1 = (Ifc+dIfc)*Vfc.case1_sym;
Pfc_case2 = (Ifc+dIfc)*Vfc.case2_sym;
Pfc.case1 = matlabFunction(Pfc_case1,'Vars',{x_normalized,Ifc_normalized,dIfc_normalized});
Pfc.case2 = matlabFunction(Pfc_case2,'Vars',{x_normalized,Ifc_normalized,dIfc_normalized});

z_fun_sym = c23*(x(4)-x(1)-x(2)-c2)/(c24*(Ifc+dIfc));
z_fun = matlabFunction(z_fun_sym,'Vars',{x_normalized,Ifc_normalized,dIfc_normalized});
% % disp(dx_case1(equil(1:4),equil(5),equil(6)))
% syms Vcm Ifc
% x = sym("x",[4 1]);
% dx = dx_case1(x,Vcm,Ifc)
% syms Vcm Ifc
% % dx = dx_case2(E,equil(4),150)
%     % x = [15042.653199799134912021585002073;
%     %   76139.944952232211555005391396959;
%     %   13655.826454600053592633027437894;
%     %   154050.0106207491901541006957663]
%     % Vcm= 236.99856301999062877802947879073
%     % Ifc= 89.527833662786594374308636991492
% 
% dx = dx_case1(x,2.316211449256888e+02,1.294426715809049e+02);
% equation = [dx(1)==0;dx(2)==0;dx(3)==0;dx(4)==0];
% 
% ax = [x];
% 
% Ans = zeros(length(ax),1);
% for ii = 1:1000
% % Answer = vpasolve(equation,ax,[0,5*patm;0,5*patm;0,1e7;0,4*patm;100,300],"Random",true)
% Answer = vpasolve(equation,ax,[1.120639379025403e+04+unifrnd(-1e4,1e4);7.167461461722918e+04+unifrnd(-7e4,7e4);1.379112974453103e+04+unifrnd(-1e4,1e4);1.442208766131414e+05+unifrnd(-1e5,1e5)])
% if isempty(Answer.x1)
%     Ans(:) = zeros(length(ax),1);
% else
% Ans(:) = double([Answer.x1;Answer.x2;Answer.x3;Answer.x4])
% end
% if (-(50*Weight_flow(Ans(3),Ans(4))-0.1)+Ans(4)/params.patm)< 0 && (-Ans(4)/params.patm+(15.27*Weight_flow(Ans(3),Ans(4))+0.6)) < 0 && Ans(1) > 1.121e+04
%     disp(Ans(:))
%     equil = Ans(:);
%     save value_equil_4d_2_2 equil
%     break
% end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Pcm = @(x,Vcm) Vcm*(Vcm-kv*x(3))/(Rcm*eata_cm);
% 
% Vfc.case1 = @(x,Ifc) Vstack1(x,Ifc,Tst,psat,patm,Tatm,Afc,tm,n);
% Vfc.case2 = @(x,Ifc) Vstack2(x,Ifc,Tst,psat,patm,Tatm,Afc,tm,n);
% 
% z_fun = @(x,Ifc) c23*(x(4)-x(1)-x(2)-c2)/(c24*Ifc);
end

function Vfc = Vstack1(x,Ifc,Tfc,psat,patm,Tatm,Afc,tm,n)
po2 = x(1);
pH2 = x(1)+x(2);


E = 1.229-0.85e-3*(Tfc-Tatm)+4.3085e-5*Tfc*(log(pH2/patm)+1/2*log(po2/patm));
v0 = 0.279-8.5e-4*(Tfc-Tatm)+4.308e-5*Tfc*(log((x(1)+x(2))/patm)+1/2*log(0.1173*(x(1)+x(2))/patm));
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
pH2 = x(1)+x(2);


E = 1.229-0.85e-3*(Tfc-Tatm)+4.3085e-5*Tfc*(log(pH2/patm)+1/2*log(po2/patm));
v0 = 0.279-8.5e-4*(Tfc-Tatm)+4.308e-5*Tfc*(log((x(1)+x(2))/patm)+1/2*log(0.1173*(x(1)+x(2))/patm));
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

%if 1e5*(po2/0.1173+psat) < 2*patm
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
