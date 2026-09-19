function [du1,K0,u_bar,x_next,u_b,ub1,Pb_bar,lambda,solve_time] = solve_snapshot(x,u_bar,Pd,params,Ref,Smax,u_b,lambda)
% One fixed OCP. This same source is installed in the two Simulink charts.
N = params.MPC.N;
du=u_bar-u_b; nc=15;
K0=zeros(3,8); X=zeros(8,N+1);
Altro.mu=repmat(4e-4,nc,N+1);
Altro.I_mu=zeros(nc,nc*(N+1));
Altro.rammda=zeros(nc,N+1);
Altro.cmax=1e-8; Altro.phi=2;
C=ones(nc,N+1);
t0=tic;
for iter=1:params.N_iter
    for k=1:N+1
        for j=1:nc
            if C(j,k)<0 && Altro.rammda(j,k)==0
                Altro.I_mu(j,nc*(k-1)+j)=0;
            else
                Altro.I_mu(j,nc*(k-1)+j)=Altro.mu(j,k);
            end
        end
    end
    [X,A,B]=forwardpass_exp(x,u_b,du,params.tsp_str,N,params,Ref);
    [du,K0,X]=backwardpass_exp(X,u_b,du,A,B,N,params.J_tsp,params.tsp_str, ...
        1,Pd,params,Altro,Ref,Smax,u_b);
    for k=1:N+1
        j=min(k,N);
        if (X(1,k)/.1173+params.psat/Ref.pres)<2
            C(:,k)=constraints_case1(X(:,k),u_b(:,j),du(:,j),Pd(k),Smax);
        else
            C(:,k)=constraints_case2(X(:,k),u_b(:,j),du(:,j),Pd(k),Smax);
        end
        % All 15 constraints are inequalities. Update on violations as well.
        Altro.rammda(:,k)=max(0,Altro.rammda(:,k)+Altro.mu(:,k).*C(:,k));
        for j=1:nc
            if C(j,k)>1e-6
                Altro.mu(j,k)=min(100,2*Altro.mu(j,k));
            end
        end
    end
end
lambda=circshift(Altro.rammda,-1,2); lambda(:,end)=lambda(:,end-1);
x_next=X(:,2); du1=du(:,1); ub1=u_b(:,1);
Pb_bar=Pb_fun(X(:,1),u_b(:,1),du(:,1));
u_b=u_b+du; u_bar=[u_b(:,2:end),u_b(:,end)];
if isfield(params,'closed_loop_feedback') && params.closed_loop_feedback==0
    K0(:)=0; % same intersample policy as the fmincon reference
end
solve_time=toc(t0);
end
