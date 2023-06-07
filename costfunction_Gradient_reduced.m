function [J,Jgrad_temp] = costfunction_Gradient_reduced(n_obs,freq,X,obs,B_r,XB_r,Sr_TL,R,H_r,G_BTQ_TL)

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
%
%  Output:
%    J:          Cost function
%    Jgradient:  Cost function gradient
%    JHessian:   Cost function Hessian
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X_red       = zeros(length(XB_r),size(X,2));
% compute reduced simulated model forecast
for i = 1:size(X,2)
    X_red(:,i)  = Sr_TL'*X(:,i);
end

% compute distances between observation and forecast
d           = zeros(size(obs,1),n_obs);
d_noisy     = zeros(size(obs,1),n_obs);
for j = n_obs:-1:1
    d(:,j)          = obs(:,j) - H_r * X_red(:,(j-1)*freq+2);
    d_noisy(:,j)    = R\(d(:,j));
    % because H is the identitiy
end
meas        = d_noisy(:);

% gradient computation - observation based
Jgrad_temp  = G_BTQ_TL'*meas;
  
% XB_r(:,1) is reduced background state (forecast) at t = 0
% X0_r = X_r(:,1) is current iterate (reduced in dim) at t = 0
% Calculate cost function - initital condition-based part
J           = 0.5*(X_red(:,1)-XB_r)'*(B_r\(X_red(:,1)-XB_r));    

% Calculate cost function - observation based part
% note indices: no observations at time t = 0!
% X and XB include time t = 0 whereas obs starts at the first time step
for i = 1:n_obs
    JOpart = d(:,i)'*d_noisy(:,i);
    J = J + 0.5 * JOpart;
end      

end