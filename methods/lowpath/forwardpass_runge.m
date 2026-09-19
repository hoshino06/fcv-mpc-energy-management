function [x_k,A,B] = forwardpass_runge(x_curr,u,du,Ist,dIst,dt,N,mode,params,Ref)
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

for ii = 1:N
    if mode == 'V'
        [~,x_phi] = runge_Kutta(@(t,x_phi) Dphi(t,Phi_u(x_phi,du(:,ii),nx),u(:,ii),Ist,dIst,params,Ref),dt,0,[X0;phi0]);
        x_k(:,ii+1) = x_phi(1:nx);
        % if X0(4)==0 && X0(9)==0
        %    x_k(4,ii+1) = 0;
        % end

        for jj=1:nx
            A_Material(:,jj) = reshape(x_phi(n*jj+1:n*(jj+1)-nu),nx,[]);
        end
        for jj=nx+1:n
            B_Material(:,jj-nx) = reshape(x_phi(n*jj+1:n*(jj+1)-nu),nx,[]);
        end

    elseif mode == 'D'
        % [~,x] = runge_Kutta(@(t,x) func(t,x,u(:,ii),du(:,ii),params,Ref),dt,0,x_k(:,ii));
        % x_k(:,ii+1) = x(end,:);
        % jac = Jac_Fun(x_k(:,ii),du(:,ii));
        % A_Material = (eye(nx)+dt*jac(1:nx,1:nx));
        % B_Material = dt*jac(1:nx,nx+1:nx+nu);
    end


    A(:,nx*(ii-1)+1:nx*ii) = A_Material;
    B(:,nu*(ii-1)+1:nu*ii) = B_Material;
    if ii == N
        break
    end
X0 = [x_k(:,ii+1);du(:,ii+1)];
end

end

function dphi = Dphi(t,phi,u,Ist,dIst,params,Ref)
% [phi(1:8),phi(9:11)] = satulation(phi(1:8),u,phi(9:11),Ref);
[phi(1:4),phi(5)] = satulation(phi(1:4),u,phi(5),Ref);
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
        dphi = case1_VE_fun(phi,u,Ist,dIst);
    else
        dphi = case2_VE_fun(phi,u,Ist,dIst);
    end
end


end

function dx = func(t,x,u,du,params,Ref)
% [x,u] = satulation(x,u,du,Ref);
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

[x,u] = satulation(x,u,du,Ref);
if (params.c11/Ref.pres-(x(1)+x(2)+params.c2/Ref.pres))>0
    dx = zeros(size(x));
else
    if params.c11/(x(1)+x(2)+params.c2/Ref.pres)>params.c19*Ref.pres||x(4)<1

        dx = Air_supply_4d_bat_case1(x,u,du);

    else
        dx = Air_supply_4d_bat_case2(x,u,du);

    end
end

end

function x_phi = Phi_u(x_phi,u,nx)
x_phi(nx+1:nx+length(u))=u;
end
