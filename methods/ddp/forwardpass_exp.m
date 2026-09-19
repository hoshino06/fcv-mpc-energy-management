function [x_k,A,B] = forwardpass_exp(x_curr,u,du,tsp_str,N,params,Ref)
                   % forwardpass_exp(x_curr,u,dt,tsp_str,N,funs,mode)
                       % forwardpass_exp     (x,u_b,du,tsp,N,D_method,params,Ref)
%UNTITLED この関数の概要をここに記述
%   詳細説明をここに記述
nx = size(x_curr,1);
nu = size(u,1);
A = zeros(nx,nx*N);
A_Material = zeros(nx);
B = zeros(nx,nu*N);
B_Material = zeros(nx,nu);
X0 = [x_curr(:,1);zeros(nu,1)];
n = length(X0);
phi0 = reshape(eye(n),[n*n,1]);
x_k=zeros(nx,N+1);
x_k(:,1) = x_curr(:,1);
% options = odeset('RelTol',1e-10,'AbsTol',1e-12);

mode = 'V'; % Vで固定
for ii = 1:N
    if mode == 'V' %変分方程式

        %%% 離散化法の切替 (params.disc_mode, 既定 0) %%%%%%%%%%%%
        %  0 = VE      : 変分方程式(次元 n+n^2)を積分して A,B を厳密に得る（本手法・売り）
        %  1 = FD中心  : 公称8状態のみ積分し、A,B は離散フローの中心差分で近似（比較用）
        %  2 = FD前進  : 同上を前進差分で近似（比較用・さらに安価）
        %  積分器(公称状態/VE とも)は params.use_rk4 で ode45 / 固定ステップ RK4 を切替。
        if isfield(params,'disc_mode'), dmode = params.disc_mode; else, dmode = 0; end
        if dmode == 0
            x_phi_end = ve_step([X0;phi0],tsp_str(ii),u(:,ii),du(:,ii),nx,params,Ref);
            x_k(:,ii+1) = x_phi_end(1:nx);
            for jj=1:nx
                A_Material(:,jj) = reshape(x_phi_end(n*jj+1:n*(jj+1)-nu),nx,[]);
            end
            for jj=nx+1:n
                B_Material(:,jj-nx) = reshape(x_phi_end(n*jj+1:n*(jj+1)-nu),nx,[]);
            end
        else
            [xe,A_Material,B_Material] = jac_step(X0(1:nx),tsp_str(ii), ...
                u(:,ii),du(:,ii),nx,nu,params,Ref,dmode==2);
            x_k(:,ii+1) = xe;
        end
        % if ii <= N_ends || ii > N-N_ends %N_endsの両端でdt=tsp_endsで動く
        % 
        %     [~,x_phi] = ode45(@(t,x_phi) Dphi(t,Phi_u(x_phi,u(:,ii),nx),funs),[0 tsp_ends],[X0;phi0]);
        %     x_k(:,ii+1) = x_phi(end,1:nx);
        % 
        %     for jj=1:nx
        %         A_Material(:,jj) = reshape(x_phi(end,n*jj+1:n*(jj+1)-nu),nx,[]);
        %     end
        %     for jj=nx+1:n
        %         B_Material(:,jj-nx) = reshape(x_phi(end,n*jj+1:n*(jj+1)-nu),nx,[]);
        %     end
        % elseif ii > N_ends || ii <= N-N_ends %予測区間の途中はdt=tsp_midで動く
        % 
        % 
        %     [~,x_phi] = ode45(@(t,x_phi) Dphi(t,Phi_u(x_phi,u(:,ii),nx),funs),[0 tsp_mid],[X0;phi0]);
        %     x_k(:,ii+1) = x_phi(end,1:nx);
        % 
        %     for jj=1:nx
        %         A_Material(:,jj) = reshape(x_phi(end,n*jj+1:n*(jj+1)-nu),nx,[]);
        %     end
        %     for jj=nx+1:n
        %         B_Material(:,jj-nx) = reshape(x_phi(end,n*jj+1:n*(jj+1)-nu),nx,[]);
        %     end
    
    % elseif mode == 'D' %差分法
    %     [~,x] = ode45(@(t,x) func(t,x,u(ii),funs),[0 dt],x_k(:,ii));
    %     x_k(:,ii+1) = x(end,:);
    %     jac = funs.Jac_Fun(x_k(:,ii),u(ii));
    %     A_Material = (eye(nx)+dt*jac(1:nx,1:nx));
    %     B_Material = dt*jac(1:nx,nx+1:nx+nu);
    end


    A(:,nx*(ii-1)+1:nx*ii) = A_Material;
    B(:,nu*(ii-1)+1:nu*ii) = B_Material;
    if ii == N
        break
    end
X0 = [x_k(:,ii+1);u(:,ii+1)];
end


% for ii = 1:N
%     if mode == 'V'
%         [~,x_phi] = ode45(@(t,x_phi) Dphi(t,Phi_u(x_phi,u(:,ii),nx),funs),[0 dt],[X0;phi0]);
%         x_k(:,ii+1) = x_phi(end,1:nx);
% 
%         for jj=1:nx
%             A_Material(:,jj) = reshape(x_phi(end,n*jj+1:n*(jj+1)-nu),nx,[]);
%         end
%         for jj=nx+1:n
%             B_Material(:,jj-nx) = reshape(x_phi(end,n*jj+1:n*(jj+1)-nu),nx,[]);
%         end
% 
%     elseif mode == 'D'
%         [~,x] = ode45(@(t,x) func(t,x,u(ii),funs),[0 dt],x_k(:,ii));
%         x_k(:,ii+1) = x(end,:);
%         jac = funs.Jac_Fun(x_k(:,ii),u(ii));
%         A_Material = (eye(nx)+dt*jac(1:nx,1:nx));
%         B_Material = dt*jac(1:nx,nx+1:nx+nu);
%     end
% 
% 
%     A(:,nx*(ii-1)+1:nx*ii) = A_Material;
%     B(:,nu*(ii-1)+1:nu*ii) = B_Material;
%     if ii == N
%         break
%     end
% X0 = [x_k(:,ii+1);u(:,ii+1)];
% end
end


% 2024年度のときのやつ
% function dphi = Dphi(t,phi,funs)
% dphi = funs.VE_Fun(phi);
% end
% 
% function dx = func(t,x,u,funs)
% dx = funs.fun(x,u);
% end
% 
% function x_phi = Phi_u(x_phi,u,nx)
% x_phi(nx+1:nx+length(u))=u;
% end

function y_end = ve_step(y0,T,u_ii,du_ii,nx,params,Ref)
% 変分方程式(次元 n+n^2)を区間 [0,T] で積分し，終端値を行ベクトルで返す。
% params.use_rk4 ~= 0 なら固定ステップ RK4，それ以外は ode45（従来と同一）。
% params.integ_h: RK4 のマイクロステップ幅[s]（既定 0.005，元 runge_Kutta.m と同じ）。
odefun = @(tt,y) Dphi(tt,Phi_u(y,du_ii,nx),u_ii,params,Ref);
use_rk4 = isfield(params,'use_rk4') && params.use_rk4 ~= 0;
if use_rk4
    if isfield(params,'integ_h'), h = params.integ_h; else, h = 0.005; end
    M  = max(1,ceil(T/h - 1e-9));
    dt = T/M;
    y  = y0(:);
    for s = 1:M
        k1 = odefun(0,y);
        k2 = odefun(0,y + 0.5*dt*k1);
        k3 = odefun(0,y + 0.5*dt*k2);
        k4 = odefun(0,y + dt*k3);
        y  = y + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
    end
    y_end = y(:).';
else
    [~,x_phi] = ode45(odefun,[0 T],y0(:));
    y_end = x_phi(end,:);
end
end

function [x_end,A_M,B_M] = jac_step(x0,T,u_ii,du_ii,nx,nu,params,Ref,fwd_diff)
% 比較用: 変分方程式を使わず、離散フロー Phi(x0,du0) の Jacobian を差分で求める。
%   Phi = 非線形8状態を区間[0,T]で RK4 積分（VE の公称状態と同一）。
%   A = dPhi/dx0, B = dPhi/ddu0 を有限差分で近似する。
%     fwd_diff=false : 中心差分（(Phi(+e)-Phi(-e))/2e）— より正確・2倍のコスト
%     fwd_diff=true  : 前進差分（(Phi(+e)-Phi(0))/e）      — 安価・1次精度
%   使う ODE ソルバ/刻みは VE と同じ（params.use_rk4, params.integ_h）。
if isfield(params,'integ_h'), h = params.integ_h; else, h = 0.005; end

x_end = flowN(x0(:),du_ii(:),T,h,u_ii,params,Ref);   % 公称（行ベクトル）
phi0  = x_end.';                                       % 列
A_M = zeros(nx,nx);
B_M = zeros(nx,nu);

sx = max(abs(x0(:)),1);   ex = 1e-6*sx;               % 状態別の摂動幅
su = max(abs(du_ii(:)),1); eu = 1e-6*su;

for j = 1:nx
    dp = zeros(nx,1); dp(j) = ex(j);
    fp = flowN(x0(:)+dp,du_ii(:),T,h,u_ii,params,Ref).';
    if fwd_diff
        A_M(:,j) = (fp - phi0) / ex(j);
    else
        fm = flowN(x0(:)-dp,du_ii(:),T,h,u_ii,params,Ref).';
        A_M(:,j) = (fp - fm) / (2*ex(j));
    end
end
for j = 1:nu
    dp = zeros(nu,1); dp(j) = eu(j);
    fp = flowN(x0(:),du_ii(:)+dp,T,h,u_ii,params,Ref).';
    if fwd_diff
        B_M(:,j) = (fp - phi0) / eu(j);
    else
        fm = flowN(x0(:),du_ii(:)-dp,T,h,u_ii,params,Ref).';
        B_M(:,j) = (fp - fm) / (2*eu(j));
    end
end
end

function xr = flowN(x0,du0,T,h,u_ii,params,Ref)
% 8状態 func を区間 [0,T] で RK4 積分し終端状態を行ベクトルで返す。
odefun = @(tt,x) func(tt,x,u_ii,du0,params,Ref);
use_rk4 = isfield(params,'use_rk4') && params.use_rk4 ~= 0;
if use_rk4
    M  = max(1,ceil(T/h - 1e-9));
    dt = T/M;
    x  = x0(:);
    for s = 1:M
        k1 = odefun(0,x);
        k2 = odefun(0,x + 0.5*dt*k1);
        k3 = odefun(0,x + 0.5*dt*k2);
        k4 = odefun(0,x + dt*k3);
        x  = x + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
    end
    xr = x(:).';
else
    [~,xx] = ode45(odefun,[0 T],x0(:));
    xr = xx(end,:);
end
end

function dphi = Dphi(t,phi,u,params,Ref)
[phi(1:8),phi(9:11)] = satulation(phi(1:8),u,phi(9:11),Ref);
% if phi(1) < 0  
%    phi(1) = 0;
% elseif phi(1) > params.Vcm_u/Ref.Vcm 
%    phi(1) = params.Vcm_u/Ref.Vcm ;
% end
% 
% if phi(2) < 0 
%     phi(2) = 0;
% elseif phi(2) > 2.2*params.Afc/Ref.Ifc   
%     phi(2) = 2.2*params.Afc/Ref.Ifc;
% end
% 
% if phi(3) <-36/Ref.Ib
%     phi(3) = -36/Ref.Ib;
% elseif phi(3) > 36/Ref.Ib
%     phi(3) = 36/Ref.Ib;
% end

if (params.c11/Ref.pres-(phi(1)+phi(2)+params.c2/Ref.pres))>0
    dphi = zeros(size(phi));
else
    if params.c11/(phi(1)+phi(2)+params.c2/Ref.pres)>params.c19*Ref.pres||phi(4)<1
        dphi = case1_VE_fun(phi,u);
    else
        dphi = case2_VE_fun(phi,u);
    end
end


end

function dx = func(t,x,u,du,params,Ref)
% backwardpass_exp.m の func と同一（元の forwardpass 版は [x,u]=satulation(...) で
% u を du で上書きする不具合があったため揃えた）。jac_step / flowN が使用。
[x,du] = satulation(x,u,du,Ref);
if (params.c11/Ref.pres-(x(1)+x(2)+params.c2/Ref.pres))>0||x(4)<1
    dx = zeros(size(x));
else
    if params.c11/(x(1)+x(2)+params.c2/Ref.pres)>params.c19*Ref.pres

        dx = Air_supply_4d_bat_case1(x,u,du);

    else
        dx = Air_supply_4d_bat_case2(x,u,du);

    end
end

end

function x_phi = Phi_u(x_phi,u,nx)
x_phi(nx+1:nx+length(u))=u;
end