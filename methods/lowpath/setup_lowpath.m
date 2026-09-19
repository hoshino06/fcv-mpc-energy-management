% Called by run_lowpath with scenario in the caller workspace.
make_function_flag = 0;
if isfield(scenario,'lowpath') && isfield(scenario.lowpath,'regenerate_functions')
    make_function_flag = scenario.lowpath.regenerate_functions;
end
addpath("model")

Ref.pres = 1.01325e5./1;
Ref.x = [Ref.pres;Ref.pres;1e5;Ref.pres]./1;
Ref.P = 1e4./1;
Ref.Vcm = 300./1;
Ref.Ifc = 2.2*280./1;
Ref.Batt = [1;1;0.1]./1;
Ref.Ib = 36./1;

[Air_supply_case1,Air_supply_case2,Pcm_fun,Pfc_fun,z_fun,params] = Air_supply_model_4d(Ref);
params.battery.parallel = scenario.plant.battery_parallel;
params.battery.series = scenario.plant.battery_series;
[Battery,Pb_fun] = Battery_model(Ref,params.battery.parallel,params.battery.series);

% params.nx = 8;
% params.nu = 3;

params.nx = 4;
params.nu = 1;


idx = [];
equil = readmatrix("equil_case1.xlsx");

params.start.Vcm = 100;
params.start.Ifc = 100;
params.start.Vcm = scenario.plant.start_Vcm_idx;
params.start.Ifc = scenario.plant.start_Ifc_idx;
params.start.x = equil(4*(params.start.Vcm-1)+1:4*params.start.Vcm,params.start.Ifc).*Ref.x;

params.start.Bat = [0.5;0;0];
params.start.Bat(1) = scenario.plant.start_SOC;
params.start.int_Ifc = 0;
params.start.u = [params.start.Vcm;params.start.Ifc;0];

params.MPC.dt = 0.05;
params.MPC.dt = scenario.dt;
params.MPC.N = 10;
params.N_iter = 200;
params.Smax = 50; %Converted to units of As by multiplying with maximum current
% params.Smax = 50;

% batteryの並列直列の数の定義
% 関数の作成には影響ない。SimuLinkのブロックの中で使われる変数
% params.battery.parallel = 15;
% params.battery.series = 12;
params.battery.parallel = 15;
params.battery.series = 12;
params.battery.parallel = scenario.plant.battery_parallel;
params.battery.series = scenario.plant.battery_series;

% Qs = 100;
% R = diag([1,1,1]);

Qs = 0.15;
R = 100;
if isfield(scenario,'lowpath')
    if isfield(scenario.lowpath,'Qs'), Qs = scenario.lowpath.Qs; end
    if isfield(scenario.lowpath,'R'),  R  = scenario.lowpath.R;  end
end

params.MPC.Qs = Qs;
params.MPC.R = R;

params.w_n = 16;
params.zeta = 1;


% symblicは全体のMPCの時と同じままで関数を作成する．
x = sym("x",[8 1]);
u = sym("u",[3 1]);
du = sym("du",[3 1]);
syms Pd Smax
Vsoc = x(5);
Vs = x(6);
Vf = x(7);
Vcm = u(1);
Ifc = u(2);
Ib = u(3);

Pb = Pb_fun(x(5:7),u(3),du(3));
Pfc_case1 = Pfc_fun.case1(x(1:4),u(2),du(2));
Pfc_case2 = Pfc_fun.case2(x(1:4),u(2),du(2));
Pd_fun1 = Pfc_case1-Pcm_fun(x(1:4),u(1),du(1))+Pb;
Pd_fun2 = Pfc_case2-Pcm_fun(x(1:4),u(1),du(1))+Pb;
Const_e_case1 = Pfc_case1-Pcm_fun(x(1:4),u(1),du(1))+Pb-Pd;
Const_e_case2 = Pfc_case2-Pcm_fun(x(1:4),u(1),du(1))+Pb-Pd;

% Const_case1 = [
%     % Const_e_case1-1;
%     % -Const_e_case1-1;
%     (params.c11/Ref.pres-(x(1)+x(2)+params.c2/Ref.pres));
%     % -Vfc_fun.case1(x(1:4),u(2)+du(2));
%     % (-Pfc_case1+Pcm_fun(x(1:4),Vcm));
%     (-(50*Weight_flow_normalized(x(3),x(4),Ref.x(3),Ref.pres)-0.1)+x(4));
%     (-x(4)+(15.27*Weight_flow_normalized(x(3),x(4),Ref.x(3),Ref.pres)+0.6));
%     (1.5-z_fun(x(1:4),u(2),du(2)));
%     (-x(1)*Ref.pres)/Ref.pres;
%     1-x(4);
%     -x(5);
%     x(5)-1;
%     (x(8)-Smax)./Smax;
%     -(u(1)+du(1));
%     ((u(1)+du(1))-params.Vcm_u/Ref.Vcm);
%     -(u(2)+du(2));
%     ((u(2)+du(2))-2.2*params.Afc/Ref.Ifc);
%     (-36/Ref.Ib-(u(3)+du(3)));
%     ((u(3)+du(3))-36/Ref.Ib);
%         % Const_e_case1/Ref.P
%     ];
Const_case1 = [
    % Const_e_case1-1;
    % -Const_e_case1-1;
    (params.c11/Ref.pres-(x(1)+x(2)+params.c2/Ref.pres));
    % -Vfc_fun.case1(x(1:4),u(2)+du(2));
    % (-Pfc_case1+Pcm_fun(x(1:4),Vcm));
    (-(50*Weight_flow_normalized(x(3),x(4),Ref.x(3),Ref.pres)-0.1)+x(4));
    (-x(4)+(15.27*Weight_flow_normalized(x(3),x(4),Ref.x(3),Ref.pres)+0.6));
    (-x(1)*Ref.pres)/Ref.pres;
    1-x(4);
    -(u(1)+du(1));
    ((u(1)+du(1))-params.Vcm_u/Ref.Vcm);
    % Const_e_case1/Ref.P
    ];
% Const_case1 = Const_case1;
Const_x_case1 = jacobian(Const_case1,x(1:4));
Const_u_case1 = jacobian(Const_case1,u(1));

% Const_case2 = [
%     % Const_e_case2-1;
%     % -Const_e_case2-1;
%     (params.c11/Ref.pres-(x(1)+x(2)+params.c2/Ref.pres));
%     % -Vfc_fun.case2(x(1:4),(u(2)+du(2)));
%     % -(Pfc_case2-Pcm_fun(x(1:4),Vcm));
%     (-(50*Weight_flow_normalized(x(3),x(4),Ref.x(3),Ref.pres)-0.1)+x(4));
%     (-x(4)+(15.27*Weight_flow_normalized(x(3),x(4),Ref.x(3),Ref.pres)+0.6));
%     (1.5-z_fun(x(1:4),u(2),du(2)));
%     (-x(1)*Ref.pres)/Ref.pres;
%     -x(5);
%     1-x(4);
%     x(5)-1;
%     (x(8)-Smax)./Smax;
%     -(u(1)+du(1));
%     ((u(1)+du(1))-params.Vcm_u/Ref.Vcm);
%     -(u(2)+du(2));
%     ((u(2)+du(2))-2.2*params.Afc/Ref.Ifc);
%     (-36/Ref.Ib-(u(3)+du(3)));
%     ((u(3)+du(3))-36/Ref.Ib);
%         % Const_e_case2/Ref.P
%     ];
Const_case2 = [
    % Const_e_case2-1;
    % -Const_e_case2-1;
    (params.c11/Ref.pres-(x(1)+x(2)+params.c2/Ref.pres));
    % -Vfc_fun.case2(x(1:4),(u(2)+du(2)));
    % -(Pfc_case2-Pcm_fun(x(1:4),Vcm));
    (-(50*Weight_flow_normalized(x(3),x(4),Ref.x(3),Ref.pres)-0.1)+x(4));
    (-x(4)+(15.27*Weight_flow_normalized(x(3),x(4),Ref.x(3),Ref.pres)+0.6));
    (-x(1)*Ref.pres)/Ref.pres;
    1-x(4);
    -(u(1)+du(1));
    ((u(1)+du(1))-params.Vcm_u/Ref.Vcm);
    % Const_e_case2/Ref.P
    ];
% Const_case2 = Const_case2;
Const_x_case2 = jacobian(Const_case2,x(1:4));
Const_u_case2 = jacobian(Const_case2,u(1));

dx_case1_sym = [Air_supply_case1(x(1:4),u(1),du(1),u(2),du(2));Battery(x(5:7),u(3),du(3));u(3)+du(3)];
dx_case2_sym = [Air_supply_case2(x(1:4),u(1),du(1),u(2),du(2));Battery(x(5:7),u(3),du(3));u(3)+du(3)];
params.nc = length(Const_case1);
params.nc_i = length(Const_case1);


%%%%%%%%%%%%%%%% IP_pendで変更したところ　%%%%%%%%%%%%%%%%%%%%%%%%
params.ns = 1;
params.MPC.mu = ones(params.nc, 1)/1e3;
%params.MPC.mu(14:15) = [1/10, 1/10];
params.N_iter_IP = 40;
params.mu_iter   = 10;
params.mu_decrease_rate = ones(params.nc,1)*5;
params.start.s = zeros(params.nc,params.MPC.N);

params.lineserch = 1;       % 1:lineserch, 0:Nonlineserch
params.D_method = 0;        % 1:differential, 0:variational

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% syms Ist dIst
Ist = sym('Ist');   
dIst = sym('dIst'); 

if make_function_flag == 1
dx_case1 = matlabFunction(dx_case1_sym,'File','functions/Air_supply_4d_bat_case1','Vars',{x,u,du});
dx_case2 = matlabFunction(dx_case2_sym,'File','functions/Air_supply_4d_bat_case2','Vars',{x,u,du});
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% lowpathのMPCに向けたfunの作成
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
dx_case1_Air = matlabFunction(Air_supply_case1(x(1:4),u(1),du(1),Ist,dIst),'File','functions/dx_case1_Air','Vars',{x,u,du,Ist,dIst});
dx_case2_Air = matlabFunction(Air_supply_case2(x(1:4),u(1),du(1),Ist,dIst),'File','functions/dx_case2_Air','Vars',{x,u,du,Ist,dIst});
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
matlabFunction(Air_supply_case1(x(1:4),u(1),du(1),u(2),du(2)),'File','functions/Air_supply_case1','Vars',{x,u,du});
matlabFunction(Air_supply_case2(x(1:4),u(1),du(1),u(2),du(2)),'File','functions/Air_supply_case2','Vars',{x,u,du});
matlabFunction(Battery(x(5:7),u(3),du(3)),'File','functions/Battery','Vars',{x(5:7),u(3),du(3)});

matlabFunction(Pd_fun1,'File','functions/Pd_fun_case1','Vars',{x,u,du});
matlabFunction(Pd_fun2,'File','functions/Pd_fun_case2','Vars',{x,u,du});
matlabFunction(Pb,'File','functions/Pb_fun','Vars',{x,u,du});
matlabFunction(z_fun(x(1:4),u(2),du(2)),'File','functions/z_fun','Vars',{x(1:4),u(2),du(2)});

matlabFunction(Const_case1,'File','functions/constraints_case1','Vars',{x,u,du,Pd,Smax})
matlabFunction(Const_x_case1,'File','functions/constraints_x_case1','Vars',{x,u,du,Pd,Smax});
matlabFunction(Const_u_case1,'File','functions/constraints_u_case1','Vars',{x,u,du,Pd,Smax});

matlabFunction(Const_case2,'File','functions/constraints_case2','Vars',{x,u,du,Pd,Smax})
matlabFunction(Const_x_case2,'File','functions/constraints_x_case2','Vars',{x,u,du,Pd,Smax});
matlabFunction(Const_u_case2,'File','functions/constraints_u_case2','Vars',{x,u,du,Pd,Smax});
addpath("functions");
counter = 0;
% [cost_case1,exfuns_case1,counter] = Cost_setting(dx_case1,params.nx,params.nu,Qs,R,Pd_fun1,Pd_fun2,counter,Ref,Pcm_fun);
% [cost_case2,exfuns_case2,counter] = Cost_setting(dx_case2,params.nx,params.nu,Qs,R,Pd_fun1,Pd_fun2,counter,Ref,Pcm_fun);

[cost_case1,exfuns_case1,counter] = Cost_setting(dx_case1_Air,params.nx,params.nu,Qs,R,Pd_fun1,Pd_fun2,counter,Ref,Pcm_fun,z_fun);
[cost_case2,exfuns_case2,counter] = Cost_setting(dx_case2_Air,params.nx,params.nu,Qs,R,Pd_fun1,Pd_fun2,counter,Ref,Pcm_fun,z_fun);

end
addpath("functions");
x_normalized = [params.start.x;params.start.Bat]./[Ref.x;Ref.Batt];
u_normalized = params.start.u./[Ref.Vcm;Ref.Ifc;Ref.Ib];
disp((x_normalized(1)/0.1173+params.psat/Ref.pres) < 2)
disp((params.start.x(1)/0.1173+params.psat) < 2*Ref.pres)
if (x_normalized(1)/0.1173+params.psat/Ref.pres) < 2
    Pd_ = matlabFunction(Pd_fun1,'Vars',{x,u,du});
    disp('case1')
else
    Pd_ = matlabFunction(Pd_fun2,'Vars',{x,u,du});
    disp('case2')
end


%% Pdの計算
% 階段状のステップ
params.Pd = Pd_(x_normalized,u_normalized,[0;0;0]);
% % open_system('make_demand.slx')
% % Data = sim("make_demand.slx");
% % close_system('make_demand.slx')
% time = 0:params.MPC.dt:20;
% Pd_time =  repmat(params.Pd,length(time),1);
% Pd_time(2/params.MPC.dt+1:end) = Pd_time(2/params.MPC.dt+1:end)+10000;
% Pd_time(2.2/params.MPC.dt+1:end) = Pd_time(2.2/params.MPC.dt+1:end)+10000;
% Data.demand = timeseries(Pd_time,time);
% Demandseries = Data.demand;
% for ii = 1:params.MPC.N
% Demandseries.Data(:,ii) = circshift(Data.demand.Data,-ii+1,1);
% end
% % Demandseries.Data = fliplr(Demandseries.Data);
% % disp(cost_case2.J_fun(x_normalized,u_normalized,params.Pd))

% Load the identical traction-demand contribution used by all B cases.
% Placeholder only; run_lowpath constructs the selected interval demand.
Data.demand = timeseries([params.Pd;params.Pd], [0;params.MPC.dt]);

time = Data.demand.Time;
Pd_base = Data.demand.Data;
% ===== MPC 用需要系列(予測区間内の需要データ) =====
Demandseries = Data.demand;
for ii = 1:params.MPC.N
    Demandseries.Data(:,ii) = circshift(Data.demand.Data, -ii+1, 1);
end


%% Pcmの初期値の計算
params.Pcm_initial = Pcm_fun(x_normalized(1:4), u_normalized(1), [0]); % 引数はPcm_funの定義に合わせる必要あり

% t=0 の Pd(0) は計算済み
params.Pd_initial = params.Pd; 

% t=0 での合計入力 P(0)
P_initial = params.Pd_initial + params.Pcm_initial;

% 積分ブロックの初期値
LPF_e_initial = [P_initial; 0];

%%
function [cost,exfuns,counter] = Cost_setting(fun,nx,nu,Qs,R,Pd_fun1,Pd_fun2,counter,Ref,Pcm_fun,z_fun_sym)
counter = counter +1;
% % % % % % % % % % % % % 
% nx = 4;
% nu = 1;
% % % % % % % % % % % % 
n = nx+nu;
phi = sym("phi",[n+n^2 1]);
u = sym("u",[nu 1]);
phi_mat = reshape(phi(n+1:end),n,n);
x = phi(1:nx);
du = phi(nx+1:nx+nu);
X = [x;du];
% syms Ist dIst
Ist = sym('Ist');   
dIst = sym('dIst'); 
% Fun = [fun(x,u,du);zeros(nu,1)];
Fun = [fun(x,u,du,Ist,dIst);zeros(nu,1)];
jac_Fun = jacobian(Fun,X);
% jac_jac_Fun = jacobian(jac_Fun,X);
dphi = [Fun;reshape(jac_Fun*phi_mat,[],1)];
exfuns.fun = fun;
exfuns.Jac_Fun = matlabFunction(jac_Fun,'File',['functions/case',num2str(counter),'_Jac_fun'],'Vars',{x,u,du,Ist,dIst});
exfuns.VE_Fun = matlabFunction(dphi,'File',['functions/case',num2str(counter),'_VE_fun'],'Vars',{phi,u,Ist,dIst});

% setting cost function
% % % % % % % % % % % % % % % % % % % % % % % % % 
% 全体のMPC
% % % % % % % % % % % % % % % % % % % % % % % % % 
% syms Pd
% du = sym('du',[nu 1]); % これ必要？
% func = ['Pd_fun_case',num2str(counter),'(x,u,du)'];
% Pd_fun = eval(func);
% % u_ref = [0;127/Ref.Ifc;0];
% % u_ref = [0;0;0];
% % u_ref = u_bar;
% J = ((Pd_fun-Pd)/Ref.P).'*Qs*((Pd_fun-Pd)/Ref.P)+du.'*R*du+1e-2*(u(2)+du(2));
% % J = (Pd_fun/Pcm_fun(x(1:4),u(1))).'*Qs*(Pd_fun/Pcm_fun(x(1:4),u(1)))+(u-u_ref).'*R*(u-u_ref)+1e-2*u(2);

% % % % % % % % % % % % % % % % % % % % % % % % % 
% LowPath MPC
% % % % % % % % % % % % % % % % % % % % % % % % % 
syms z
du = sym('du',[nu 1]);
z_fun = z_fun_sym(x,Ist,dIst);


J = (z_fun-z)*Qs*(z_fun-z) + du*R*du; 

lx = jacobian(J,x);
lxx = jacobian(lx,x);
lu = jacobian(J,du);
luu = jacobian(lu,du);
lux = jacobian(lu,x);
%  全体の時のMPC
% J_fun = matlabFunction(J,'File',['functions/case',num2str(counter),'_J_fun'],'Vars',{x,u,du,Pd});
% lx_fun = matlabFunction(lx,'File',['functions/case',num2str(counter),'_lx_fun'],'Vars',{x,u,du,Pd});
% lxx_fun = matlabFunction(lxx,'File',['functions/case',num2str(counter),'_lxx_fun'],'Vars',{x,u,du,Pd});
% lu_fun = matlabFunction(lu,'File',['functions/case',num2str(counter),'_lu_fun'],'Vars',{x,u,du,Pd});
% luu_fun = matlabFunction(luu,'File',['functions/case',num2str(counter),'_luu_fun'],'Vars',{x,u,du,Pd});
% lux_fun = matlabFunction(lux,'File',['functions/case',num2str(counter),'_lux_fun'],'Vars',{x,u,du,Pd});

% lowpathの時のMPC
J_fun = matlabFunction(J,'File',['functions/case',num2str(counter),'_J_fun'],'Vars',{x,du,Ist,dIst,z});
lx_fun = matlabFunction(lx,'File',['functions/case',num2str(counter),'_lx_fun'],'Vars',{x,du,Ist,dIst,z});
lxx_fun = matlabFunction(lxx,'File',['functions/case',num2str(counter),'_lxx_fun'],'Vars',{x,du,Ist,dIst,z});
lu_fun = matlabFunction(lu,'File',['functions/case',num2str(counter),'_lu_fun'],'Vars',{x,du,Ist,dIst,z});
luu_fun = matlabFunction(luu,'File',['functions/case',num2str(counter),'_luu_fun'],'Vars',{x,du,Ist,dIst,z});
lux_fun = matlabFunction(lux,'File',['functions/case',num2str(counter),'_lux_fun'],'Vars',{x,du,Ist,dIst,z});


funs = who('-regexp','fun$');
for i = 1:size(funs,1)
    cost.(funs{i}) = eval(funs{i});
end
end

% bus_setting;
% bus_setting_old;
