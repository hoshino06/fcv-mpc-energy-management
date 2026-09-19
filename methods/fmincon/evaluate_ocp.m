function [J,c,X] = evaluate_ocp(U,x,Uref,pd,params,Ref,Smax)
% Shared definition for same-state solver comparisons. Inputs are normalized.
N=params.MPC.N; U=reshape(U,3,N); steps=params.tsp_str;
X=zeros(8,N+1); X(:,1)=x; c=zeros(params.nc*(N+1),1); J=0;
for k=1:N+1
    j=min(k,N); u=Uref(:,j); du=U(:,j)-u;
    if (x(1)/.1173+params.psat/Ref.pres)<2
        ell=case1_J_fun(x,u,du,pd(k));
        ci=constraints_case1(x,u,du,pd(k),Smax);
    else
        ell=case2_J_fun(x,u,du,pd(k));
        ci=constraints_case2(x,u,du,pd(k),Smax);
    end
    if k<=N, w=params.J_tsp(k); else, w=params.terminal_weight; end
    J=J+w*ell;
    c((k-1)*params.nc+(1:params.nc))=ci;
    if k<=N
        n=max(1,ceil(steps(k)/params.integ_h)); h=steps(k)/n;
        for q=1:n
            a=fcv_dyn(x,U(:,k),params,Ref);
            b=fcv_dyn(x+h*a/2,U(:,k),params,Ref);
            d=fcv_dyn(x+h*b/2,U(:,k),params,Ref);
            e=fcv_dyn(x+h*d,U(:,k),params,Ref);
            x=x+h*(a+2*b+2*d+e)/6;
        end
        X(:,k+1)=x;
    end
end
end
