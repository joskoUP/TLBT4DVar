function [J,Jgradient,JHessian] = costfunction_fixed_H(tstep,h,XB,X,obs,B,R,H,freq,force)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Purpose: Calculate cost function 
%  
%  Input:
%    h:            Time step for numerical scheme
%    freq:         Frequency of observations
%    tstep:        Number of time steps to perform
%    XB:           Background solution (forecast)
%    X:            Current state; forward trajectory as rows 
%    observations: Observation values as rows
%    force:        parameter for Lorenz model ffunc
%
%  Output:
%    J:          Cost function
%    Jgradient:  Cost function gradient
%    JHessian:   Cost function Hessian
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% XB(:,1) is background state (forecast) at t = 0
% X0 = X(:,1) is current iterate at t = 0
% TLBT only works for nonzero B, so assume this here
J       = 0.5*(X(:,1)-XB)'*(B\(X(:,1)-XB));    

n_obs   = ceil(tstep/freq); 

% Calculate cost function
% note indices: no observations at time t = 0!
% X and XB include time t = 0 whereas obs starts at the first time step
for i = 1:n_obs
    JOpart  = (obs(:,i)-H*X(:,(i-1)*freq+2))'*(R\(obs(:,i)-H*X(:,(i-1)*freq+2)));
    J       = J + 0.5 * JOpart;
end

% Calculate the adjoint
Xadj = zeros(size(B,1),1);
GradXadj = zeros(size(B));

% this does not work for positive semidefinite matrices
% L = chol((H'*(R\H)),'lower');
% this works for positive semidefinite matrices
[L,D] = ldl((H'*(R\H)),'lower');
L = L*sqrt(D);

for j = n_obs:-1:1
    
    % Calculate M_{j,j-1}(X_{j-1})
    M = JacobianRK4(h*freq,X(:,(j-1)*freq+1),force);   
       
    Xadj = M'*(Xadj + H'*(R\(obs(:,j) - H*X(:,(j-1)*freq+2))));
    GradXadj = [GradXadj; L'] * M;
 
end                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      

Jgradient =(B\(X(:,1)-XB)) - Xadj;
JHessian = (B\eye(size(B))) + GradXadj'*GradXadj;

end