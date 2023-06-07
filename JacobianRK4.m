function M = JacobianRK4(h,state,force)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Find Jacobian of the solver (RK4 method) to 
%  obtain matrix M_{n+1,n}
%
%  2008  M. A. Freitag
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Variables:
%
%
%  Input:
%
%    h:          Time step for numerical scheme
%    state:      state variable, which includes
%                [Qinput, Pinput] Initial fields (INPUT)
%   force:       parameter for Lorenz model ffunc
%
%  Output:
%
%    M:          Jacobian M_{n+1,n} 
%                from one time step to the next one
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X0      = state;

s       = size(state,1);

% Implement RK method
 
X1      = X0 + 0.5*h*ffunc(X0,force);
 
X2      = X0 + 0.5*h*ffunc(X1,force);
    
X3      = X0 + h*ffunc(X2,force);

fpz1    = [fgradfunc(X0)];
fpz2    = [fgradfunc(X1)];
fpz3    = [fgradfunc(X2)];
fpz4    = [fgradfunc(X3)];

ifpz1   = eye(s) + h*fpz1/2;
ifpz2   = eye(s) + h*fpz2*ifpz1/2;
ifpz3   = eye(s) + h*fpz3*ifpz2;

M       = eye(s) + h/6 * (fpz1 + 2*fpz2*ifpz1 + 2*fpz3*ifpz2 + fpz4*ifpz3);

end


