function [statetraj,hall] = rk4(tstep,h,stateinitial,force)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Apply 4th order Runge-Kutta to the Lorenz model
%
%  2008  M. A. Freitag
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Variables:
%
%  Input:
%
%    h:            Time step for numerical scheme
%    tstep:        Number of time steps to perform
%    stateinitial: Initial fields (INPUT)
%    force:        parameter for Lorenz model ffunc
%
%  Output:
%
%    statetraj:    Trajectories of evolved fields (OUTPUT)
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N           = length (stateinitial);

X           = zeros(N,tstep);

X(:,1)      = stateinitial;

% Implement RK method

hall        = h;

for i = 1:tstep
        
    X0          = X(:,i);    
 
    X1          = X0 + 0.5*h*ffunc(X0,force);
        
    X2          = X0 + 0.5*h*ffunc(X1,force);
    
    X3          = X0 + h*ffunc(X2,force);
            
    X(:,i+1)    = X0 + h*(ffunc(X0,force)+2*ffunc(X1,force)+ ...
                    2*ffunc(X2,force)+ffunc(X3,force))/6; 
  
end

% output

statetraj   = X;

end