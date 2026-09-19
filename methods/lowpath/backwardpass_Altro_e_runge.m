function [u_opt,x_opt,K0,l_check] = backwardpass_Altro_e_runge(x_bar,u_bar,du,Ist,dIst,A,B,N,dt,linesearch,Pd,params,Altro,Ref,Smax,z,u_ref)

if linesearch == 1
    N_alpha = 5;
    alpha = 0.5.^(0:N_alpha-1);
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
% l_check
l_check = zeros(1,N+1);

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
            lx(:,ii) = case1_lx_fun(x_bar(:,ii),du(:,ii),Ist,dIst,z).';%+DBdx(x_bar(:,ii),Smax);
            lxx(:,nx*(ii-1)+1:nx*ii) = case1_lxx_fun(x_bar(:,ii),du(:,ii),Ist,dIst,z);%+DBdxx(x_bar(:,ii),Smax) ;
            lu(:,ii) = case1_lu_fun(x_bar(:,ii),du(:,ii),Ist,dIst,z).';
            luu(:,nu*(ii-1)+1:nu*ii) = case1_luu_fun(x_bar(:,ii),du(:,ii),Ist,dIst,z);
            Jbs(ii) = case1_J_fun(x_bar(:,ii),du(:,ii),Ist,dIst,z);
            lux(:,nx*(ii-1)+1:nx*ii) = case1_lux_fun(x_bar(:,ii),du(:,ii),Ist,dIst,z);
            
            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
            % lの値を見るために追加する
            % disp([z, z_fun(x_bar(:,ii),Ist,dIst)])
            l_check(:,ii) = case1_J_fun(x_bar(:,ii),0,Ist,dIst,z);
            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
        else
            % lx(:,ii) = case2_lx_fun(x_bar(:,ii),u_bar(:,ii)).';
            % lxx(:,nx*(ii-1)+1:nx*ii) = case2_lxx_fun(x_bar(:,ii),u_bar(:,ii)) ;
            % lu(:,ii) = case2_lu_fun(x_bar(:,ii),u_bar(:,ii)).';
            % luu(:,nu*(ii-1)+1:nu*ii) = case2_luu_fun(x_bar(:,ii),u_bar(:,ii));
            lx(:,ii) = case2_lx_fun(x_bar(:,ii),du(:,ii),Ist,dIst,z).';%+DBdx(x_bar(:,ii),Smax);
            lxx(:,nx*(ii-1)+1:nx*ii) = case2_lxx_fun(x_bar(:,ii),du(:,ii),Ist,dIst,z);%+DBdxx(x_bar(:,ii),Smax) ;
            lu(:,ii) = case2_lu_fun(x_bar(:,ii),du(:,ii),Ist,dIst,z).';
            luu(:,nu*(ii-1)+1:nu*ii) = case2_luu_fun(x_bar(:,ii),du(:,ii),Ist,dIst,z);
            lux(:,nx*(ii-1)+1:nx*ii) = case2_lux_fun(x_bar(:,ii),du(:,ii),Ist,dIst,z);
            Jbs(ii) = case2_J_fun(x_bar(:,ii),du(:,ii),Ist,dIst,z);

            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
            % lの値を見るために追加する
            l_check(:,ii) = case2_J_fun(x_bar(:,ii),0,Ist,dIst,z);
            % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
        end

        Jb = sum(Jbs);

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
        lx(:,ii) = case1_lx_fun(x_bar(:,ii),du(:,N),Ist,dIst,z).';%+DBdx(x_bar(:,ii),Smax);
        lxx(:,nx*(ii-1)+1:nx*ii) = case1_lxx_fun(x_bar(:,ii),du(:,N),Ist,dIst,z);%+DBdxx(x_bar(:,ii),Smax) ;
        % lx(:,ii) = zeros(nx,1);
        % lxx(:,nx*(ii-1)+1:nx*ii) = zeros(nx);

        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
        % lの値を見るために追加する
        l_check(:,ii) = case1_J_fun(x_bar(:,ii),0,Ist,dIst,z);
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

    else
        lx(:,ii) = case2_lx_fun(x_bar(:,ii),du(:,N),Ist,dIst,z).';%+DBdx(x_bar(:,ii),Smax);
        lxx(:,nx*(ii-1)+1:nx*ii) = case2_lxx_fun(x_bar(:,ii),du(:,N),Ist,dIst,z);%+DBdxx(x_bar(:,ii),Smax) ;
        % lx(:,ii) = zeros(nx,1);
        % lxx(:,nx*(ii-1)+1:nx*ii) = zeros(nx);

        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
        % lの値を見るために追加する
        l_check(:,ii) = case2_J_fun(x_bar(:,ii),0,Ist,dIst,z);
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    end
    if (x_bar(1,ii)/0.1173+params.psat/Ref.pres) < 2
        c(:,ii) = constraints_case1(x_bar(:,ii),u_bar(:,N),du(:,N),Pd(ii-1),Smax);
        c_x(:,nx*(ii-1)+1:nx*ii)  = constraints_x_case1(x_bar(:,ii),u_bar(:,N),du(:,N),Pd(N),Smax);
        c_u(:,nu*(ii-1)+1:nu*ii)  = constraints_u_case1(x_bar(:,ii),u_bar(:,N),du(:,N),Pd(N),Smax);
    else
        c(:,ii) = constraints_case2(x_bar(:,ii),u_bar(:,N),du(:,N),Pd(ii-1),Smax);
        c_x(:,nx*(ii-1)+1:nx*ii)  = constraints_x_case2(x_bar(:,ii),u_bar(:,N),du(:,N),Pd(N),Smax);
        c_u(:,nu*(ii-1)+1:nu*ii)  = constraints_u_case2(x_bar(:,ii),u_bar(:,N),du(:,N),Pd(N),Smax);
    end

    p(:,end) = lx(:,end)+c_x(:,end-nx+1:end) .'*(rammda(:,end)+I_mu(:,end-nc+1:end)*c(:,end));
    P(:,end-nx+1:end) = lxx(:,end-nx+1:end) + c_x(:,end-nx+1:end).'*I_mu(:,nc*(ii-1)+1:nc*ii)*c_x(:,end-nx+1:end);
    % if any(isinf(c))
    %     disp(c)
    % end


    for ii=N:-1:1

        Qxx = lxx(:,nx*(ii-1)+1:nx*ii)+A(:,nx*(ii-1)+1:nx*(ii)).'*P(:,nx*ii+1:nx*(ii+1))*A(:,nx*(ii-1)+1:nx*(ii))+c_x(:,nx*(ii-1)+1:nx*ii) .'*I_mu(:,nc*(ii-1)+1:nc*ii)*c_x(:,nx*(ii-1)+1:nx*ii) ;
        Quu = luu(:,nu*(ii-1)+1:nu*ii)+B(:,nu*(ii-1)+1:nu*ii).'*P(:,nx*ii+1:nx*(ii+1))*B(:,nu*(ii-1)+1:nu*ii)+c_u(:,nu*(ii-1)+1:nu*ii) .'*I_mu(:,nc*(ii-1)+1:nc*ii)*c_u(:,nu*(ii-1)+1:nu*ii) ;
        % Quxにluxの値入れてないけどいいんか？
        % Qux =  B(:,nu*(ii-1)+1:nu*ii).'*P(:,nx*(ii+1)-nx+1:nx*(ii+1))*(A(:,nx*(ii)-nx+1:nx*(ii))) + c_u(:,nu*(ii-1)+1:nu*ii) .'*I_mu(:,nc*(ii-1)+1:nc*ii)*c_x(:,nx*(ii-1)+1:nx*ii) ;
        Qux = lux(:,nx*(ii-1)+1:nx*ii) + B(:,nu*(ii-1)+1:nu*ii).'*P(:,nx*(ii+1)-nx+1:nx*(ii+1))*(A(:,nx*(ii)-nx+1:nx*(ii))) + c_u(:,nu*(ii-1)+1:nu*ii) .'*I_mu(:,nc*(ii-1)+1:nc*ii)*c_x(:,nx*(ii-1)+1:nx*ii) ;
        Qxu = Qux.';
        Qx = lx(:,ii)+A(:,nx*(ii-1)+1:nx*(ii)).'*p(:,ii+1) + c_x(:,nx*(ii-1)+1:nx*ii) .'*(rammda(:,ii)+I_mu(:,nc*(ii-1)+1:nc*ii)*c(:,ii));
        Qu = lu(:,ii)+B(:,nu*(ii-1)+1:nu*ii).'*p(:,ii+1) + c_u(:,nu*(ii-1)+1:nu*ii).'*(rammda(:,ii)+I_mu(:,nc*(ii-1)+1:nc*ii)*c(:,ii));
        if det(Quu) ==0 && rho ==0
            rho = 1e-8;
        end
        K(:,:,ii) = -(Quu+ eye(nu)*rho)\Qux;
        d(:,ii) = -(Quu+ eye(nu)*rho)\Qu;

        P(:,nx*ii-nx+1:nx*ii) = Qxx+K(:,:,ii).'*Quu*K(:,:,ii)+K(:,:,ii).'*Qux+Qxu*K(:,:,ii);
        p(:,ii) = Qx +K(:,:,ii).'*Quu*d(:,ii)+K(:,:,ii).'*Qu+Qxu*d(:,ii);

    end

    xk = x_bar;
    % for ii =1:N
    % Js(ii) = J_fun(xk(:,ii),u_a(:,ii));
    % end
    % disp(sum(Js))

    % linesearch
    for jj =1:N_alpha
        % disp('-------------')
        for ii = 1:N
            u_a(:,ii) = du(:,ii)+K(:,:,ii)*(xk(:,ii)-x_bar(:,ii))+alpha(jj)*d(:,ii);
            [~,u_a(:,ii)] = satulation(x_bar(:,ii),u_bar(:,ii),u_a(:,ii),Ref);
            % disp(u_bar(:,ii)+u_a(:,ii))
            [~,x] = runge_Kutta(@(t,x) func(t,x,u_bar(:,ii),u_a(:,ii),Ist,dIst,params,Ref),dt,0,xk(:,ii));
            xk(:,ii+1) = x;

            if (xk(1,ii)/0.1173+params.psat/Ref.pres) < 2
                if any(xk(:,ii)<=0)
                Js(ii) = 1e8;
                else
                Js(ii) = case1_J_fun(xk(:,ii),u_a(:,ii),Ist,dIst,z);
                end
            else
                if any(xk(:,ii)<=0)
                Js(ii) = 1e8;
                else
                Js(ii) = case2_J_fun(xk(:,ii),u_a(:,ii),Ist,dIst,z);
                end
            end
        end
        x_cand(:,:,jj) = xk;
        u_cand(:,:,jj) = u_a;
        % if (xk(1,ii)/0.1173+params.psat/Ref.pres) < 2
        %     J_cand(jj) = sum(Js);%+case1_J_fun(xk(:,N+1),zeros(nu,1),zeros(nu,1),Pd(end));
        % else
            J_cand(jj) = sum(Js);%+case2_J_fun(xk(:,N+1),zeros(nu,1),zeros(nu,1),Pd(end));
        % end
    end
    [J,idx_opt] = min(J_cand);
    % disp(J)
    % disp('----------------')
    x_opt = x_cand(:,:,idx_opt);

    % if flag == 1
    %         for ii = 1:N
    %         u_a(:,ii) = u_bar(:,ii)+K(:,:,ii)*(xk(:,ii)-x_bar(:,ii))+d(:,ii);
    %         end
    %         u_opt = u_a;
    % else
    u_opt = u_cand(:,1:N,idx_opt);
    % end

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


function dx = func(t,x,u,du,Ist,dIst,params,Ref)
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

        % dx = Air_supply_4d_bat_case1(x,u,du);
        dx = dx_case1_Air(x,u,du,Ist,dIst);

    else
        % dx = Air_supply_4d_bat_case2(x,u,du);
        dx = dx_case2_Air(x,u,du,Ist,dIst);

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
