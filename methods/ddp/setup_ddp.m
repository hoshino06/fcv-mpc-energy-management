% Preserve the scenario supplied by the experiment runner across workspace reset.
setup_scenario = evalin('base','BM_scn');
clearvars -except setup_scenario
make_function_flag = 0;
if isfield(setup_scenario,'cost') && isfield(setup_scenario.cost,'regenerate_functions')
    make_function_flag = setup_scenario.cost.regenerate_functions;
end
addpath("model")

Ref.pres = 1.01325e5./1;
Ref.x = [Ref.pres;Ref.pres;1e5;Ref.pres]./1;
Ref.P = 1e4./1;
Ref.Vcm = 300./1;
Ref.Ifc = 2.2*280./1;
Ref.Batt = [1;1;0.1]./1;
Ref.Ib = 36./1;
Max_Ib = setup_scenario.limits.Ib_abs_max;
% Ref.Ib = 400./1;
% Max_Ib = 400;

[Air_supply_case1,Air_supply_case2,Pcm_fun,Pfc_fun,z_fun,params] = Air_supply_model_4d(Ref);
params.battery.parallel = setup_scenario.plant.battery_parallel;
params.battery.series = setup_scenario.plant.battery_series;
[Battery,Pb_fun] = Battery_model(Ref,params.battery.parallel,params.battery.series);

params.nx = 8;
params.nu = 3;


idx = [];
equil = readmatrix("equil_case1.xlsx");

params.start.Vcm = 100;
params.start.Ifc = 100;
params.start.Vcm = setup_scenario.plant.start_Vcm_idx;
params.start.Ifc = setup_scenario.plant.start_Ifc_idx;
params.start.x = equil(4*(params.start.Vcm-1)+1:4*params.start.Vcm,params.start.Ifc).*Ref.x;

params.start.Bat = [0.5;0;0];
params.start.Bat(1) = setup_scenario.plant.start_SOC;
params.start.int_Ifc = 0;
params.start.u = [params.start.Vcm;params.start.Ifc;0];

params.MPC.dt = 0.05; % control period [s]
params.MPC.N = 5;     % non-uniform horizon steps
params.N_iter = 200;

% --- MPC 内部積分器の選択（計算時間削減 #1）----------------------------------
% use_rk4 = 0 : 従来の可変ステップ ode45（基準・数値検証用）
% use_rk4 = 1 : 固定ステップ RK4（実時間実装パス）。forwardpass_exp / backwardpass_exp が参照。
% integ_h     : RK4 のマイクロステップ幅[s]。各予測区間の分割数 = ceil(tsp_str_i / integ_h)。
% ※ Simulink チャートに params が構造体パラメータとして渡るため、文字列フィールドは
%   使わず数値フラグにすること（char フィールドは runtime parameter に写像できない）。
params.use_rk4 = 0;      % 既定は基準の ode45。実時間計測時に 1 へ。
params.integ_h = 0.025;  % RK4 刻み[s]。sweep_integ_h.m の結果: 0.025 で対 ode45 相対差 ~2e-5、約6倍高速

% 離散化法（forwardpass の A,B の求め方）。本手法は 0（変分方程式）。1/2 は比較用。
%  0 = VE     : 変分方程式を積分（次元 n+n^2）。非一様・大きい Δt でも A,B が正確（本研究の売り）
%  1 = FD中心 : 公称8状態のみ積分し、離散フローの中心差分で A,B を近似
%  2 = FD前進 : 同上を前進差分で（安価・1次精度）
params.disc_mode = 0;
params.terminal_weight = 1; % fixed terminal weight for both horizon grids
params.closed_loop_feedback = 0; % hold the first optimized input for 50 ms

% backwardpass の line search: 0 = 従来の総当たり(alpha 5 点全評価) / 1 = バックトラッキング
% alpha を大きい順に評価し、十分減少条件 (実測減少/予測減少 の比 z ∈ [c1,c2]) を満たした時点で打ち切る。
% Howell et al. の ALTRO forward pass 相当。SquareRoot 等は不使用。配列形状は総当たりと同一。
% 計測: 同一解を返しつつ line search を約1.5倍高速化（N_iter が大きいほど効く）。
params.ls_backtrack = 1;
params.ls_c1        = 1e-4;   % 十分減少比の下限
params.ls_c2        = 10;     % 同上限
params.Smax = 72/Ref.Ib; %Converted to units of As by multiplying with maximum current
% params.Smax = 50;

% batteryの並列直列数　Simulinkのブロックのみに対応する．関数作成には対応していない．
params.battery.parallel = setup_scenario.plant.battery_parallel;
params.battery.series = setup_scenario.plant.battery_series;

Qs = 100;
R = diag([1,1,0.01]);
if isfield(setup_scenario,'cost')
    if isfield(setup_scenario.cost,'Qs'), Qs=setup_scenario.cost.Qs; end
    if isfield(setup_scenario.cost,'R'), R=setup_scenario.cost.R; end
end

params.MPC.Qs = Qs;
params.MPC.R = R;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% tsp_strの作成　（不均一ホライゾン）

t_end = 0.5; % prediction horizon [s]
% 最初のステップを dt_first にするための係数 alpha を見つける 
% horizon_param.mで計算
alpha = horizon_param(params.MPC.N,t_end,params.MPC.dt);
% alpha = 3.3126; % N = 5のとき

sum_exp = 0;
tsp_str = zeros(1,params.MPC.N);

for i = 1:params.MPC.N
    sum_exp = sum_exp + exp(i/alpha);
end

for i = 1:params.MPC.N
    tsp_str(i) = exp(i/alpha) / sum_exp * t_end;
end

% 予測区間におけるステップ幅tsp_strから，予測区間の時間軸tsp_str1を計算
tsp_str1 = tsp_str;
for i = 1:size(tsp_str1,2)-1
    tsp_str1(i+1) = tsp_str1(i) + tsp_str1(i+1);
end
% 実時間のサンプル時間tspの間に入る， 実際の入力の数
% for i = 1:params.MPC.N
%     if tsp_str1(i)<tsp
%         s = i+1;
%     end
% end
J_tsp = tsp_str/0.05; % もともとのホライゾン内のステップ幅 0.05s で割る．


% simulinkのILQRブロックに送るため構造体paramに格納
params.tsp_str = tsp_str;
params.J_tsp = J_tsp;


% tsp_strを固定周期にして，元のILQRと比較
 % tsp_str = 0.05*ones(1,params.MPC.N);
 % J_tsp = tsp_str/0.05;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 





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

Const_case1 = [
    % Const_e_case1-1;
    % -Const_e_case1-1;
    (params.c11/Ref.pres-(x(1)+x(2)+params.c2/Ref.pres));
    % -Vfc_fun.case1(x(1:4),u(2)+du(2));
    % (-Pfc_case1+Pcm_fun(x(1:4),Vcm));
    (-(50*Weight_flow_normalized(x(3),x(4),Ref.x(3),Ref.pres)-0.1)+x(4));
    (-x(4)+(15.27*Weight_flow_normalized(x(3),x(4),Ref.x(3),Ref.pres)+0.6));
    (1.5-z_fun(x(1:4),u(2),du(2)));
    (-x(1)*Ref.pres)/Ref.pres;
    1-x(4);
    -x(5);
    x(5)-1;
    (x(8)-Smax)./Smax;
    -(u(1)+du(1));
    ((u(1)+du(1))-params.Vcm_u/Ref.Vcm);
    -(u(2)+du(2));
    ((u(2)+du(2))-2.2*params.Afc/Ref.Ifc);
    % (-36/Ref.Ib-(u(3)+du(3)));
    % ((u(3)+du(3))-36/Ref.Ib);
    (-Max_Ib/Ref.Ib-(u(3)+du(3)));
    ((u(3)+du(3))-Max_Ib/Ref.Ib);
        % Const_e_case1/Ref.P
    ];
% Const_case1 = Const_case1;
Const_x_case1 = jacobian(Const_case1,x);
Const_u_case1 = jacobian(Const_case1,u);

Const_case2 = [
    % Const_e_case2-1;
    % -Const_e_case2-1;
    (params.c11/Ref.pres-(x(1)+x(2)+params.c2/Ref.pres));
    % -Vfc_fun.case2(x(1:4),(u(2)+du(2)));
    % -(Pfc_case2-Pcm_fun(x(1:4),Vcm));
    (-(50*Weight_flow_normalized(x(3),x(4),Ref.x(3),Ref.pres)-0.1)+x(4));
    (-x(4)+(15.27*Weight_flow_normalized(x(3),x(4),Ref.x(3),Ref.pres)+0.6));
    (1.5-z_fun(x(1:4),u(2),du(2)));
    (-x(1)*Ref.pres)/Ref.pres;
    -x(5);
    1-x(4);
    x(5)-1;
    (x(8)-Smax)./Smax;
    -(u(1)+du(1));
    ((u(1)+du(1))-params.Vcm_u/Ref.Vcm);
    -(u(2)+du(2));
    ((u(2)+du(2))-2.2*params.Afc/Ref.Ifc);
    % (-36/Ref.Ib-(u(3)+du(3)));
    % ((u(3)+du(3))-36/Ref.Ib);
    (-Max_Ib/Ref.Ib-(u(3)+du(3)));
    ((u(3)+du(3))-Max_Ib/Ref.Ib);
        % Const_e_case2/Ref.P
    ];
% Const_case2 = Const_case2;
Const_x_case2 = jacobian(Const_case2,x);
Const_u_case2 = jacobian(Const_case2,u);

dx_case1_sym = [Air_supply_case1(x(1:4),u(1),du(1),u(2),du(2));Battery(x(5:7),u(3),du(3));u(3)+du(3)];
dx_case2_sym = [Air_supply_case2(x(1:4),u(1),du(1),u(2),du(2));Battery(x(5:7),u(3),du(3));u(3)+du(3)];
params.nc = length(Const_case1);
params.nc_i = length(Const_case1);



if make_function_flag == 1
dx_case1 = matlabFunction(dx_case1_sym,'File','functions/Air_supply_4d_bat_case1','Vars',{x,u,du});
dx_case2 = matlabFunction(dx_case2_sym,'File','functions/Air_supply_4d_bat_case2','Vars',{x,u,du});

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
[cost_case1,exfuns_case1,counter] = Cost_setting(dx_case1,params.nx,params.nu,Qs,R,Pd_fun1,Pd_fun2,counter,Ref,Pcm_fun);
[cost_case2,exfuns_case2,counter] = Cost_setting(dx_case2,params.nx,params.nu,Qs,R,Pd_fun1,Pd_fun2,counter,Ref,Pcm_fun);

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
% params.Pd = Pd_(x_normalized,u_normalized,[0;0;0]);
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

%% 不均一ホライゾンPd
% params.MPC.N = 5; % いったん再定義
params.Pd = Pd_(x_normalized,u_normalized,[0;0;0]);
time = 0:params.MPC.dt:20;
Pd_base =  repmat(params.Pd,length(time),1);

% 階段状のステップにする
% Pd_base(2/params.MPC.dt+1:end) = Pd_base(2/params.MPC.dt+1:end)+10000;
% Pd_base(2.2/params.MPC.dt+1:end) = Pd_base(2.2/params.MPC.dt+1:end)+10000;
% simulation_time = 5;

% ドライビングサイクルから求めた電力需要変化
% ===== PFt付きドライビングサイクルの読み込み =====
file_path = fullfile('..','..','studies','03_udds','data',"都市ダイナモメーター運転スケジュール_PFt付き.txt");
cycle_data = readtable(file_path,'VariableNamingRule','preserve');
time_cycle_all = cycle_data{:,1};      % Time [s]（1秒刻み）
PFt_cycle_all  = cycle_data{:,end};    % PFt [W]

% ===== 必要な区間のみ抽出 =====567, 582
% idx_start = 575;
% idx_end   = idx_start+6;
% % idx_end   = 580;

idx_start = 567;
idx_end   = 582;
simulation_time = idx_end-idx_start; % simulinkのシミュレーション時間

time_cycle = time_cycle_all(idx_start:idx_end);
PFt_cycle  = PFt_cycle_all(idx_start:idx_end);

% ===== MPCのdtに合わせて再サンプリング =====
time_cycle_dt = (time_cycle(1):params.MPC.dt:time_cycle(end)).'; % dt刻みの時間軸を作成
Pd_time = interp1(time_cycle, PFt_cycle, ...
                  time_cycle_dt, 'linear'); % 1次ホールド
Pd_time(isnan(Pd_time)) = 0;       % 念のため NaN はゼロ埋め
Data.demand = timeseries( params.Pd + Pd_time/5 , time_cycle_dt - time_cycle_dt(1)); 
                                            % ^^^^
                                            % ここが小さすぎてSmaxが上がりきらなかったex)1/50じゃだめ．

time = Data.demand.Time;
Pd_base = Data.demand.Data;

% --- 不均一ホライゾン用の多列データの作成 ---
% 各行(時刻t)において、将来の [t+tsp_str1(1), t+tsp_str1(2), ..., t+tsp_str1(N)] の値を取得する
Pd_future = zeros(length(time), params.MPC.N+1);

eps_time = 1e-4; % 微小時間を追加することで内挿のエラーを防ぐ
lookahead = [0 tsp_str1];
for ii = 1:params.MPC.N + 1
    % 今から何秒後の値が必要か
    lookahead_time = lookahead(ii) + eps_time; 
    
    % interp1 を使用して、ZOH(0次ホールド)を再現 ('previous')
    % time_base に対して、(現在の時刻 + 将来へのオフセット) でサンプリングする
    Pd_future(:, ii) = interp1(time, Pd_base, time + lookahead_time, 'previous', 'extrap');
end

% Simulinkに渡すためのtimeseries作成
Data.demand = timeseries(Pd_future, time);
Demandseries = Data.demand;



%% Pcmの初期値の計算
params.Pcm_initial = Pcm_fun(x_normalized(1:4), u_normalized(1), [0]); % 引数はPcm_funの定義に合わせる必要あり

% t=0 の Pd(0) は計算済み
params.Pd_initial = params.Pd; 

% t=0 での合計入力 P(0)
P_initial = params.Pd_initial + params.Pcm_initial;

% 積分ブロックの初期値
LPF_e_initial = [P_initial; 0];

%%
function [cost,exfuns,counter] = Cost_setting(fun,nx,nu,Qs,R,Pd_fun1,Pd_fun2,counter,Ref,Pcm_fun)
counter = counter +1;
n = nx+nu;
phi = sym("phi",[n+n^2 1]);
u = sym("u",[nu 1]);
phi_mat = reshape(phi(n+1:end),n,n);
x = phi(1:nx);
du = phi(nx+1:nx+nu);
X = [x;du];
Fun = [fun(x,u,du);zeros(nu,1)];
jac_Fun = jacobian(Fun,X);
% jac_jac_Fun = jacobian(jac_Fun,X);
dphi = [Fun;reshape(jac_Fun*phi_mat,[],1)];
exfuns.fun = fun;
exfuns.Jac_Fun = matlabFunction(jac_Fun,'File',['functions/case',num2str(counter),'_Jac_fun'],'Vars',{x,u,du});
exfuns.VE_Fun = matlabFunction(dphi,'File',['functions/case',num2str(counter),'_VE_fun'],'Vars',{phi,u});

% setting cost function
syms Pd
du = sym('du',[nu 1]);
func = ['Pd_fun_case',num2str(counter),'(x,u,du)'];
Pd_fun = eval(func);
% u_ref = [0;127/Ref.Ifc;0];
% u_ref = [0;0;0];
% u_ref = u_bar;
J = ((Pd_fun-Pd)/Ref.P).'*Qs*((Pd_fun-Pd)/Ref.P)+du.'*R*du+1e-2*(u(2)+du(2));
% J = (Pd_fun/Pcm_fun(x(1:4),u(1))).'*Qs*(Pd_fun/Pcm_fun(x(1:4),u(1)))+(u-u_ref).'*R*(u-u_ref)+1e-2*u(2);
lx = jacobian(J,x);
lxx = jacobian(lx,x);
lu = jacobian(J,du);
luu = jacobian(lu,du);
lux = jacobian(lu,x);

J_fun = matlabFunction(J,'File',['functions/case',num2str(counter),'_J_fun'],'Vars',{x,u,du,Pd});
lx_fun = matlabFunction(lx,'File',['functions/case',num2str(counter),'_lx_fun'],'Vars',{x,u,du,Pd});
lxx_fun = matlabFunction(lxx,'File',['functions/case',num2str(counter),'_lxx_fun'],'Vars',{x,u,du,Pd});
lu_fun = matlabFunction(lu,'File',['functions/case',num2str(counter),'_lu_fun'],'Vars',{x,u,du,Pd});
luu_fun = matlabFunction(luu,'File',['functions/case',num2str(counter),'_luu_fun'],'Vars',{x,u,du,Pd});
lux_fun = matlabFunction(lux,'File',['functions/case',num2str(counter),'_lux_fun'],'Vars',{x,u,du,Pd});


funs = who('-regexp','fun$');
for i = 1:size(funs,1)
    cost.(funs{i}) = eval(funs{i});
end
end

% bus_setting;
% bus_setting_old;

%% --- 初期状態の表示 ---
fprintf('\n==================================================\n');
fprintf('   Simulation Initial Settings (Steady State)\n');
fprintf('==================================================\n');

% 1. インデックスと入力値
fprintf('[Inputs]\n');
fprintf('  Vcm Index: %d (Value: %.1f V)\n', params.start.Vcm, params.start.Vcm); % インデックス=電圧の場合
fprintf('  Ifc Index: %d (Value: %.1f A)\n', params.start.Ifc, params.start.Ifc);
fprintf('  Bat SOC:   %.2f\n', params.start.Bat(1));

% 2. 電力関係 (W)
fprintf('\n[Power Flows]\n');
fprintf('  Demand Power (Pd):  %.2f W\n', params.Pd);
fprintf('  Compressor Power:   %.2f W\n', params.Pcm_initial);
fprintf('  Total FC Output:    %.2f W\n', params.Pd + params.Pcm_initial);

% 3. 物理状態 (内部変数)
z_init = z_fun(x_normalized(1:4), u_normalized(2), 0); % 酸素過剰率の計算
P_manifold = params.start.x(4) / 1000; % kPaに変換
omega_rad = params.start.x(3);              % rad/s?

fprintf('\n[System Status]\n');
fprintf('  Oxygen Excess Ratio (z): %.3f\n', z_init);
fprintf('  Manifold Pressure:       %.2f kPa (abs)\n', P_manifold);
fprintf('  Compressor Speed (w_cm):  %.2f rad/s ? \n', omega_rad);

% 4. 適用モデルの確認
if (x_normalized(1)/0.1173 + params.psat/Ref.pres) < 2
    fprintf('  Applied Model:           Case 1 (Subsonic)\n');
else
    fprintf('  Applied Model:           Case 2 (Choked)\n');
end
fprintf('==================================================\n\n');
