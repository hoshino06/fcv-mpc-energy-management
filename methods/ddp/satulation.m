function [x,du] = satulation(x,u,du,Ref)
%UNTITLED2 この関数の概要をここに記述
%   詳細説明をここに記述
if u(1)+du(1) < 0  
   % du(1) = 0-u(1);
   du(1) = 0-u(1);
elseif u(1)+du(1) > 1 
   % du(1) = 1-u(1) ;
   du(1) = 1-u(1) ;
end

if u(2)+du(2) < 0 
    % du(2) = 0-u(2);
    du(2) = 0-u(2);
elseif u(2)+du(2) > 1  
    % du(2) = 1-u(2);
    du(2) = 1-u(2);
end

if u(3)+du(3) <-1
    % du(3) = -1-u(3);
    du(3) = -1-u(3);
elseif u(3)+du(3) > 1
    % du(3) = 1-u(3);
    du(3) = 1-u(3);
end


% Battery polarization voltages and net discharged charge can be negative.
for ii=1:min(4,length(x))
    if x(ii) <=0 || isnan(x(ii))
        x(ii) = 1e-12;
    end
end
end
