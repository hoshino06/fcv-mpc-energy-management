function [u_opt,K0,x_opt] = backwardpass_exp(x_bar,u_bar,du,A,B,N,J_tsp,tsp_str,linesearch,Pd,params,Altro,Ref,Smax,u_ref)

if linesearch == 1
    N_alpha = 6;
    alpha = [0.5.^(0:4),0];
else
    N_alpha = 1;
    alpha = 1;
end


nx = size(x_bar,1);
nu = size(u_bar,1);
lx = zeros(nx,N+1);
lu = zeros(nu,N);
lxx = zeros(nx,(N+1)*nx);
luu = zeros(nu,N*nu);
lux = zeros(nu,N*nx);
p = zeros(nx,N+1);
d = zeros(nu,N);
K = zeros(nu,nx,N);
P = zeros(nx,(N+1)*nx);
rho = 0;

nc = params.nc;
nc_i= params.nc_i;
mu = Altro.mu;
I_mu = Altro.I_mu;
rammda = Altro.rammda;
c = zeros(nc,N+1);
c_x = zeros(nc,(N+1)*nx);
c_u = zeros(nc,(N+1)*nu);

x_cand = zeros(nx,N+1,N_alpha);
u_cand = zeros(nu,N,N_alpha);
J_cand = NaN(1,N_alpha);
Js = zeros(1,N);
Jbs = zeros(1,N);
u_a = u_bar;


for k = 1:30
    for ii = 1:N
        %saturation
        [x_bar(:,ii),du(:,ii)] = satulation(x_bar(:,ii),u_bar(:,ii),du(:,ii),Ref);
        % if u_bar(1,ii) < 0
        %    u_bar(1,ii) = 0;
        % elseif u_bar(1,ii) > params.Vcm_u/Ref.Vcm
        %    u_bar(1,ii) = params.Vcm_u/Ref.Vcm ;
        % end
        %
        % if u_bar(2,ii) < 0
        %     u_bar(2,ii) = 0;
        % elseif u_bar(2,ii) > 2.2*params.Afc/Ref.Ifc
        %     u_bar(2,ii) = 2.2*params.Afc/Ref.Ifc;
        % end
        %
        % if u_bar(3,ii) <-36/Ref.Ib
        %     u_bar(3,ii) = -36/Ref.Ib;
        % elseif u_bar(3,ii) > 36/Ref.Ib
        %     u_bar(3,ii) = 36/Ref.Ib;
        % end

        if (x_bar(1,ii)/0.1173+params.psat/Ref.pres) < 2
            % lx(:,ii) = case1_lx_fun(x_bar(:,ii),u_bar(:,ii)).';
            % lxx(:,nx*(ii-1)+1:nx*ii) = case1_lxx_fun(x_bar(:,ii),u_bar(:,ii)) ;
            % lu(:,ii) = case1_lu_fun(x_bar(:,ii),u_bar(:,ii)).';
            % luu(:,nu*(ii-1)+1:nu*ii) = case1_luu_fun(x_bar(:,ii),u_bar(:,ii));
            lx(:,ii) = case1_lx_fun(x_bar(:,ii),u_bar(:,ii),du(:,ii),Pd(ii)).' * J_tsp(ii);%+DBdx(x_bar(:,ii),Smax);
            lxx(:,nx*(ii-1)+1:nx*ii) = case1_lxx_fun(x_bar(:,ii),u_bar(:,ii),du(:,ii),Pd(ii)) * J_tsp(ii);%+DBdxx(x_bar(:,ii),Smax) ;
            lu(:,ii) = case1_lu_fun(x_bar(:,ii),u_bar(:,ii),du(:,ii),Pd(ii)).' * J_tsp(ii);
            luu(:,nu*(ii-1)+1:nu*ii) = case1_luu_fun(x_bar(:,ii),u_bar(:,ii),du(:,ii),Pd(ii)) * J_tsp(ii);
            Jbs(ii) = case1_J_fun(x_bar(:,ii),u_bar(:,ii),du(:,ii),Pd(ii));
            lux(:,nx*(ii-1)+1:nx*ii) = case1_lux_fun(x_bar(:,ii),u_bar(:,ii),du(:,ii),Pd(ii)) * J_tsp(ii);
        else
            % lx(:,ii) = case2_lx_fun(x_bar(:,ii),u_bar(:,ii)).';
            % lxx(:,nx*(ii-1)+1:nx*ii) = case2_lxx_fun(x_bar(:,ii),u_bar(:,ii)) ;
            % lu(:,ii) = case2_lu_fun(x_bar(:,ii),u_bar(:,ii)).';
            % luu(:,nu*(ii-1)+1:nu*ii) = case2_luu_fun(x_bar(:,ii),u_bar(:,ii));
            lx(:,ii) = case2_lx_fun(x_bar(:,ii),u_bar(:,ii),du(:,ii),Pd(ii)).' * J_tsp(ii);%+DBdx(x_bar(:,ii),Smax);
            lxx(:,nx*(ii-1)+1:nx*ii) = case2_lxx_fun(x_bar(:,ii),u_bar(:,ii),du(:,ii),Pd(ii)) * J_tsp(ii);%+DBdxx(x_bar(:,ii),Smax) ;
            lu(:,ii) = case2_lu_fun(x_bar(:,ii),u_bar(:,ii),du(:,ii),Pd(ii)).' * J_tsp(ii);
            luu(:,nu*(ii-1)+1:nu*ii) = case2_luu_fun(x_bar(:,ii),u_bar(:,ii),du(:,ii),Pd(ii)) * J_tsp(ii);
            lux(:,nx*(ii-1)+1:nx*ii) = case2_lux_fun(x_bar(:,ii),u_bar(:,ii),du(:,ii),Pd(ii)) * J_tsp(ii);
            Jbs(ii) = case2_J_fun(x_bar(:,ii),u_bar(:,ii),du(:,ii),Pd(ii));
        end

        Jb = sum(Jbs.*J_tsp);

        if (x_bar(1,ii)/0.1173+params.psat/Ref.pres) < 2

            c(:,ii) = constraints_case1(x_bar(:,ii),u_bar(:,ii),du(:,ii),Pd(ii),Smax);
            c_x(:,nx*(ii-1)+1:nx*ii)  = constraints_x_case1(x_bar(:,ii),u_bar(:,ii),du(:,ii),Pd(ii),Smax);
            c_u(:,nu*(ii-1)+1:nu*ii)  = constraints_u_case1(x_bar(:,ii),u_bar(:,ii),du(:,ii),Pd(ii),Smax);
        else
            c(:,ii) = constraints_case2(x_bar(:,ii),u_bar(:,ii),du(:,ii),Pd(ii),Smax);
            c_x(:,nx*(ii-1)+1:nx*ii)  = constraints_x_case2(x_bar(:,ii),u_bar(:,ii),du(:,ii),Pd(ii),Smax);
            c_u(:,nu*(ii-1)+1:nu*ii)  = constraints_u_case2(x_bar(:,ii),u_bar(:,ii),du(:,ii),Pd(ii),Smax);

        end

        % if c(:,ii) < 0 && rammda(:,ii) == zeros(nc,1)
        %     I_mu(:,nc*(ii-1)+1:nc*ii) = zeros(nc);
        % else
        %     I_mu(:,nc*(ii-1)+1:nc*ii) = diag(mu);
        % end

    end

    % final condition
    ii=N+1;
    [x_bar(:,ii),~] = satulation(x_bar(:,ii),zeros(nu,1),zeros(nu,1),Ref);
    if (x_bar(1,ii)/0.1173+params.psat/Ref.pres) < 2
        % lx(:,ii) = case1_lx_fun(x_bar(:,ii),u_bar(:,N));
        % lxx(:,nx*(ii-1)+1:nx*ii) = case1_lxx_fun(x_bar(:,ii),u_bar(:,N)) ;
        lx(:,ii) = case1_lx_fun(x_bar(:,ii),u_bar(:,N),du(:,N),Pd(ii)).' * params.terminal_weight;%+DBdx(x_bar(:,ii),Smax);
        lxx(:,nx*(ii-1)+1:nx*ii) = case1_lxx_fun(x_bar(:,ii),u_bar(:,N),du(:,N),Pd(ii)) * params.terminal_weight;%+DBdxx(x_bar(:,ii),Smax) ;
        % lx(:,ii) = zeros(nx,1);
        % lxx(:,nx*(ii-1)+1:nx*ii) = zeros(nx);
    else
        lx(:,ii) = case2_lx_fun(x_bar(:,ii),u_bar(:,N),du(:,N),Pd(ii)).' * params.terminal_weight;%+DBdx(x_bar(:,ii),Smax);
        lxx(:,nx*(ii-1)+1:nx*ii) = case2_lxx_fun(x_bar(:,ii),u_bar(:,N),du(:,N),Pd(ii)) * params.terminal_weight;%+DBdxx(x_bar(:,ii),Smax) ;
        % lx(:,ii) = zeros(nx,1);
        % lxx(:,nx*(ii-1)+1:nx*ii) = zeros(nx);
    end
    if (x_bar(1,ii)/0.1173+params.psat/Ref.pres) < 2
        c(:,ii) = constraints_case1(x_bar(:,ii),u_bar(:,N),du(:,N),Pd(ii),Smax);
        c_x(:,nx*(ii-1)+1:nx*ii)  = constraints_x_case1(x_bar(:,ii),u_bar(:,N),du(:,N),Pd(ii),Smax);
        c_u(:,nu*(ii-1)+1:nu*ii)  = constraints_u_case1(x_bar(:,ii),u_bar(:,N),du(:,N),Pd(ii),Smax);
    else
        c(:,ii) = constraints_case2(x_bar(:,ii),u_bar(:,N),du(:,N),Pd(ii),Smax);
        c_x(:,nx*(ii-1)+1:nx*ii)  = constraints_x_case2(x_bar(:,ii),u_bar(:,N),du(:,N),Pd(ii),Smax);
        c_u(:,nu*(ii-1)+1:nu*ii)  = constraints_u_case2(x_bar(:,ii),u_bar(:,N),du(:,N),Pd(ii),Smax);
    end

    Jb = Jb + merit_terminal(x_bar(:,end),u_bar(:,end),du(:,end),Pd(end),params,Ref,params.terminal_weight);
    for im=1:N+1
        Jb=Jb+rammda(:,im).'*c(:,im)+.5*c(:,im).'*I_mu(:,nc*(im-1)+(1:nc))*c(:,im);
    end
    p(:,end) = lx(:,end)+c_x(:,end-nx+1:end) .'*(rammda(:,end)+I_mu(:,end-nc+1:end)*c(:,end));
    P(:,end-nx+1:end) = lxx(:,end-nx+1:end) + c_x(:,end-nx+1:end).'*I_mu(:,nc*(ii-1)+1:nc*ii)*c_x(:,end-nx+1:end);
    % if any(isinf(c))
    %     disp(c)
    % end


    % Terminal terms also depend on the last decision input. Include that
    % direct dependence, not only the path through the terminal state.
    if (x_bar(1,end)/.1173+params.psat/Ref.pres)<2
        tu=case1_lu_fun(x_bar(:,end),u_bar(:,end),du(:,end),Pd(end)).'*params.terminal_weight;
        tuu=case1_luu_fun(x_bar(:,end),u_bar(:,end),du(:,end),Pd(end))*params.terminal_weight;
        tux=case1_lux_fun(x_bar(:,end),u_bar(:,end),du(:,end),Pd(end))*params.terminal_weight;
    else
        tu=case2_lu_fun(x_bar(:,end),u_bar(:,end),du(:,end),Pd(end)).'*params.terminal_weight;
        tuu=case2_luu_fun(x_bar(:,end),u_bar(:,end),du(:,end),Pd(end))*params.terminal_weight;
        tux=case2_lux_fun(x_bar(:,end),u_bar(:,end),du(:,end),Pd(end))*params.terminal_weight;
    end
    ctu=c_u(:,end-nu+1:end); ctx=c_x(:,end-nx+1:end); W=I_mu(:,end-nc+1:end);
    tu=tu+ctu.'*(rammda(:,end)+W*c(:,end));
    tuu=tuu+ctu.'*W*ctu; tux=tux+ctu.'*W*ctx;
    dV1 = 0; dV2 = 0;
    for ii=N:-1:1

        Qxx = lxx(:,nx*(ii-1)+1:nx*ii)+A(:,nx*(ii-1)+1:nx*(ii)).'*P(:,nx*ii+1:nx*(ii+1))*A(:,nx*(ii-1)+1:nx*(ii))+c_x(:,nx*(ii-1)+1:nx*ii) .'*I_mu(:,nc*(ii-1)+1:nc*ii)*c_x(:,nx*(ii-1)+1:nx*ii) ;
        Quu = luu(:,nu*(ii-1)+1:nu*ii)+B(:,nu*(ii-1)+1:nu*ii).'*P(:,nx*ii+1:nx*(ii+1))*B(:,nu*(ii-1)+1:nu*ii)+c_u(:,nu*(ii-1)+1:nu*ii) .'*I_mu(:,nc*(ii-1)+1:nc*ii)*c_u(:,nu*(ii-1)+1:nu*ii) ;
        Qux = lux(:,nx*(ii-1)+1:nx*ii) + B(:,nu*(ii-1)+1:nu*ii).'*P(:,nx*(ii+1)-nx+1:nx*(ii+1))*(A(:,nx*(ii)-nx+1:nx*(ii))) + c_u(:,nu*(ii-1)+1:nu*ii) .'*I_mu(:,nc*(ii-1)+1:nc*ii)*c_x(:,nx*(ii-1)+1:nx*ii) ;
        if ii==N
            At=A(:,nx*(ii-1)+(1:nx)); Bt=B(:,nu*(ii-1)+(1:nu));
            Qux=Qux+tux*At;
            Quu=Quu+tuu+Bt.'*tux.'+tux*Bt;
        end
        Qxu = Qux.';
        Qx = lx(:,ii)+A(:,nx*(ii-1)+1:nx*(ii)).'*p(:,ii+1) + c_x(:,nx*(ii-1)+1:nx*ii) .'*(rammda(:,ii)+I_mu(:,nc*(ii-1)+1:nc*ii)*c(:,ii));
        Qu = lu(:,ii)+B(:,nu*(ii-1)+1:nu*ii).'*p(:,ii+1) + c_u(:,nu*(ii-1)+1:nu*ii).'*(rammda(:,ii)+I_mu(:,nc*(ii-1)+1:nc*ii)*c(:,ii));
        if ii==N, Qu=Qu+tu; end
        if det(Quu) ==0 && rho ==0
            rho = 1e-8;
        end
        K(:,:,ii) = -(Quu+ eye(nu)*rho)\Qux;
        d(:,ii) = -(Quu+ eye(nu)*rho)\Qu;

        dV1 = dV1 + d(:,ii).'*Qu;
        dV2 = dV2 + 0.5*(d(:,ii).'*Quu*d(:,ii));

        P(:,nx*ii-nx+1:nx*ii) = Qxx+K(:,:,ii).'*Quu*K(:,:,ii)+K(:,:,ii).'*Qux+Qxu*K(:,:,ii);
        p(:,ii) = Qx +K(:,:,ii).'*Quu*d(:,ii)+K(:,:,ii).'*Qu+Qxu*d(:,ii);

    end

    xk = x_bar;
    % for ii =1:N
    % Js(ii) = J_fun(xk(:,ii),u_a(:,ii));
    % end
    % disp(sum(Js))

    %  lineserch
    % ls_bt = 0: 従来の総当たり（alpha を全点評価して最小コストを採用）
    % ls_bt = 1: バックトラッキング（alpha=1 から順に評価し、十分減少条件を満たした時点で打ち切り）
    %            ALTRO の forward pass 相当。SquareRoot 等は不使用。配列形状は総当たりと同一。
    if isfield(params,'ls_backtrack'), ls_bt = params.ls_backtrack ~= 0; else, ls_bt = false; end
    if isfield(params,'ls_c1'), ls_c1 = params.ls_c1; else, ls_c1 = 1e-4; end
    if isfield(params,'ls_c2'), ls_c2 = params.ls_c2; else, ls_c2 = 10;   end

    for jj =1:N_alpha
        % A line-search rollout is only a hypothetical DDP candidate. Keep
        % generated model functions out of their known real-valued domain.
        % The nominal trajectory and the real plant remain fail-fast.
        % disp('-------------')
        for ii = 1:N
            u_a(:,ii) = du(:,ii)+K(:,:,ii)*(xk(:,ii)-x_bar(:,ii))+alpha(jj)*d(:,ii);
            [~,u_a(:,ii)] = satulation(x_bar(:,ii),u_bar(:,ii),u_a(:,ii),Ref);
            % disp(u_bar(:,ii)+u_a(:,ii))
            xk(:,ii+1) = plant_step(xk(:,ii),tsp_str(ii),u_bar(:,ii),u_a(:,ii),params,Ref);

            if (xk(1,ii)/0.1173+params.psat/Ref.pres) < 2
                if any(xk(1:4,ii)<=0) || any(~isfinite(xk(:,ii)))
                Js(ii) = 1e8;
                else
                Js(ii) = case1_J_fun(xk(:,ii),u_bar(:,ii),u_a(:,ii),Pd(ii));
                end
            else
                if any(xk(1:4,ii)<=0) || any(~isfinite(xk(:,ii)))
                Js(ii) = 1e8;
                else
                Js(ii) = case2_J_fun(xk(:,ii),u_bar(:,ii),u_a(:,ii),Pd(ii));
                end
            end
        end
        x_cand(:,:,jj) = xk;
        u_cand(:,:,jj) = u_a;
        % if (xk(1,ii)/0.1173+params.psat/Ref.pres) < 2
        %     J_cand(jj) = sum(Js);%+case1_J_fun(xk(:,N+1),zeros(nu,1),zeros(nu,1),Pd(end));
        % else
            candidate_valid = isreal(xk) && all(isfinite(xk(:))) && ...
                all(reshape(xk(1:4,:)>0,[],1));
            if candidate_valid
            J_cand(jj) = sum(Js.*J_tsp)+merit_terminal(xk(:,end),u_bar(:,end),u_a(:,end),Pd(end),params,Ref,params.terminal_weight);
            for im=1:N+1
                iu=min(im,N);
                if (xk(1,im)/.1173+params.psat/Ref.pres)<2
                    cm=constraints_case1(xk(:,im),u_bar(:,iu),u_a(:,iu),Pd(im),Smax);
                else
                    cm=constraints_case2(xk(:,im),u_bar(:,iu),u_a(:,iu),Pd(im),Smax);
                end
                J_cand(jj)=J_cand(jj)+rammda(:,im).'*cm+.5*cm.'*I_mu(:,nc*(im-1)+(1:nc))*cm;
            end
            if ~isreal(J_cand(jj)) || ~isfinite(J_cand(jj))
                candidate_valid = false;
            end
            else
                J_cand(jj)=Inf;
            end
        % end

        if ~candidate_valid
            J_cand(jj) = Inf;
            x_cand(:,:,jj) = x_bar;
            u_cand(:,:,jj) = du;
        end

        if ls_bt
            % Howell et al. の受理条件: 実測減少 / 予測減少 の比 z が [c1,c2] に入れば採用
            expdec = -(alpha(jj)*dV1 + alpha(jj)^2*dV2);   % 予測コスト減少量(>0 を期待)
            zrat   = (Jb - J_cand(jj)) / expdec;
            ok = (J_cand(jj) < 1e8) && ...
                 ( (expdec > 0 && zrat > ls_c1 && zrat < ls_c2) || ...
                   (expdec <= 0 && J_cand(jj) < Jb) );
            if ok
                break
            end
        end
    end

    [J,idx_opt] = min(J_cand);   % 未評価の alpha は NaN のまま。min は NaN を無視する
    % disp(J)
    x_opt = x_cand(:,:,idx_opt);
    u_opt = u_cand(:,1:N,idx_opt);

    K0 = K(:,:,1);
    if k ==1
        rho = rho + 1e-8;
    elseif rho <= 1
        rho = rho*2;
    end
    if J < 1e8 
        break
    end
end
end


function x_end = plant_step(x0,T,u_ii,du_ii,params,Ref)
% 8状態プラントを区間 [0,T] で積分して終端状態(列ベクトル)を返す。
% params.use_rk4 ~= 0 で固定ステップ RK4，それ以外は従来通り ode45。
odefun = @(tt,x) func(tt,x,u_ii,du_ii,params,Ref);
use_rk4 = isfield(params,'use_rk4') && params.use_rk4 ~= 0;
if use_rk4
    if isfield(params,'integ_h'), h = params.integ_h; else, h = 0.005; end
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
    x_end = x(:);
else
    [~,x] = ode45(odefun,[0 T],x0(:));
    x_end = x(end,:).';
end
end

function dx = func(t,x,u,du,params,Ref)
% if params.c11/params.patm-(x(1)+x(2)+params.c2/params.patm) > 0
%  x(1) = params.c11/params.patm-(x(2)+params.c2/params.patm);
% end
% dx = Air_supply_4d_bat(x,u);
% if u(1) < 0
%    u(1) = 150;
% elseif u(1) > params.Vcm_u/Ref.Vcm
%    u(1) = params.Vcm_u/Ref.Vcm ;
% end
%
% if u(2) < 0
%     u(2) = 0;
% elseif u(2) > 2.2*params.Afc/Ref.Ifc
%     u(2) = 2.2*params.Afc/Ref.Ifc;
% end
%
% if u(3) <-36/Ref.Ib
%     u(3) = -36/Ref.Ib;
% elseif u(3) > 36/Ref.Ib
%     u(3) = 36/Ref.Ib;
% end

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


% function dBdx = DBdx(x,Smax)
%     delta = 1e-2; %small numel
%     zs =  -x(8)+Smax;
%     dzdx = -1;
%     dBdx = zeros(size(x));
%     for i = 1:length(zs) %拘束条件の数（＝バリア関数の数）だけ繰り返す
%         dBdz = barrier_dz(zs(i),delta);
%         dBdx(8) = dBdx(8) + dBdz*dzdx(i,:).'; %dBdzにdzduをかけて足していく
%     end
%     dBdx = dBdx *1e-8;
%     % dBdx =0;
%     disp([dBdx(8),x(8)])
% end
%
% function dBdxx = DBdxx(x,Smax)
%     delta = 1e-2; %small numel
%     zs =  -x(8)+Smax;
%     dzdx = -1;
%     dBdxx = zeros(length(x));
%
%     for i = 1:length(zs)
%         dBdzz = barrier_hes_z(zs(i),delta);
%         dBdxx = dBdxx + dzdx(i,:).'*dBdzz*dzdx(i,:);
%     end
%     dBdxx = dBdxx *1e-8;
% end
% %バリア関数の微分値（B(z)のz微分）
% function value = barrier_dz(z,delta)
%     if z > delta
%         value = -(1/z);
%     else
%         value = (z-2*delta)/delta;
%     end
% end
%
% function value = barrier_hes_z(z,delta)
%     if z > delta
%         value = 1/(z^2);
%     else
%         value = 1/delta;
%     end
% end

function J=merit_terminal(x,u,du,pd,params,Ref,w)
if (x(1)/.1173+params.psat/Ref.pres)<2
    J=w*case1_J_fun(x,u,du,pd);
else
    J=w*case2_J_fun(x,u,du,pd);
end
end
