function dx = Air_supply_4d(x,u)
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


%x = [Vsoc;Vs;Vf]
% dx = @(x,Ib) [-x(1)/(Rsd*Cb)-Ib/Cb;
%     -x(2)/(Rs(x(1)*100)*Cs(x(1)*100))+Ib/Cs(x(1)*100);
%     -x(3)/(Rf(x(1)*100)*Cf(x(1)*100))+Ib/Cf(x(1)*100)];
% load value_equil_4d.mat
% 
% E(1) = (Kapper-c5)*equil(1)+c2-c6;
% E(2) = equil(1)-E(1);
% E(3) = equil(2);
% E(4) = equil(3);


if c11/(x(1)+x(2)+c2) > c19
dx = [
    c1*(-x(1)-x(2)+x(4)-c2)-(c3*x(1)*psi_case1(x,c2,c11,c12,c17,c18,c19,c20))/(c4*x(1)+c5*x(2)+c6)-c7*u(2);
    c8*(-x(1)-x(2)+x(4)-c2)-(c3*x(2)*psi_case1(x,c2,c11,c12,c17,c18,c19,c20))/(c4*x(1)+c5*x(2)+c6);
    -c9*x(3)-(c10./x(3))*((x(4)./c11).^c12-1)*Weight_flow(x(3),x(4))+c13*u(1);
    c14*(1+c15*((x(4)./c11)^c12-1))*(Weight_flow(x(3),x(4))-c16*(-x(1)-x(2)+x(4)-c2))
    ];
else

dx =  [
    c1*(-x(1)-x(2)+x(4)-c2)-(c3*x(1)*psi_case2(x,c2,c11,c12,c17,c18,c19,c20))/(c4*x(1)+c5*x(2)+c6)-c7*u(2);
    c8*(-x(1)-x(2)+x(4)-c2)-(c3*x(2)*psi_case2(x,c2,c11,c12,c17,c18,c19,c20))/(c4*x(1)+c5*x(2)+c6);
    -c9*x(3)-(c10./x(3))*((x(4)./c11).^c12-1)*Weight_flow(x(3),x(4))+c13*u(1);
    c14*(1+c15*((x(4)./c11)^c12-1))*(Weight_flow(x(3),x(4))-c16*(-x(1)-x(2)+x(4)-c2))];
end




% % disp(dx_case1(equil(1:4),equil(5),equil(6)))
% x = sym("x",[4 1]);
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
% 
% Pcm = @(x,Vcm) Vcm*(Vcm-kv*x(3))/(Rcm*eata_cm);
% 
% Vfc.case1 = @(x,Ifc) Vstack1(x,Ifc,Tst,psat,patm,Tatm,Afc,tm,n);
% Vfc.case2 = @(x,Ifc) Vstack2(x,Ifc,Tst,psat,patm,Tatm,Afc,tm,n);
% z_fun = @(x,Ifc) c23*(x(4)-x(1)-x(2)-c2)/(c24*Ifc);
end

function psi = psi_case1(x,c2,c11,c12,c17,c18,c19,c20)
psi = c17*(x(1)+x(2)+c2)*(c11/(x(1)+x(2)+c2))^c18*sqrt(1-(c11/(x(1)+x(2)+c2))^c12);
end

function psi = psi_case2(x,c2,c11,c12,c17,c18,c19,c20)
psi = c20*(x(1)+x(2)+c2);
end
