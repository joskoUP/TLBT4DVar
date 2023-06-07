function [JHess,B_r,H_r,XB_r,Sr_TL,G_BTQ_TL] = Simple_TLBT_Quantities_lorenz(r,N,n_obs,endtime,B,R,M,H,XB)

% if rank is chosen too big, take full model
if r>N
    r = N;
end

%% set given prior = Reachability Gramian
Gamma_pr    = B;
L_pr        = chol(B);
A           = M*Gamma_pr+Gamma_pr*M';
if sum(real(eig(A))>0) > 0
    disp('The prior covariance L_pr is not prior-compatible.')
end

G = H'*(R\H);

%% compute time-limited Obs Gramian
% for unstable A, Gramian computation directly via summation
step    = endtime/n_obs;
Mat     = expm(M*step);
Q_TL = 0;
for k = 0:n_obs-1
Q_TL = Q_TL + ((Mat')^k)*G*Mat^k;
end
% compute a square root factorization of Q_TL
% floating point computation errors induce complex zeros
[L,D] = ldl(Q_TL,'lower');
L_Q_TL = real(L)*real(sqrt(D));

%% time-limited balancing with Q_TL 
[V,S,W]     = svd(L_Q_TL'*L_pr);

S           = S(1:r,1:r);
delQ_TL     = diag(S);
Siginvsqrt  = diag(1./sqrt(delQ_TL));
Sr_TL       = (Siginvsqrt*V(:,1:r)'*L_Q_TL')';
Tr_TL       = L_pr*W(:,1:r)*Siginvsqrt; % balancing transformation
M_r         = Sr_TL'*M*Tr_TL;
H_r         = H*Tr_TL;
B_r         = Sr_TL'*B*Sr_TL;
XB_r        = Sr_TL'*XB;

% Balancing with Q_TL - generate G_BT,H_BT
G_BTQ_TL        = zeros(n_obs*(size(H,1)),N+1);
G_BTQo_TL       = zeros(n_obs*(size(H,1)),N+1);
iter            = expm(M_r*step);
temp            = H_r;
for i = 1:n_obs
    temp                                            = temp*iter;
    temp_proj                                       = temp*Sr_TL';
    G_BTQ_TL((i-1)*(size(H,1))+1:i*(size(H,1)),:)   = temp_proj;
    G_BTQo_TL((i-1)*(size(H,1))+1:i*(size(H,1)),:)  = sqrt(R)\temp_proj;
end
L_prinv           = L_pr\eye(N+1);

% Balancing with Q_TL - compute posterior covariance and mean
R_posinv        = qr([G_BTQo_TL; L_prinv],0);
R_posinv        = triu(R_posinv(1:(N+1),:)); % Pull out upper triangular factor
JHess           = R_posinv' * R_posinv;

end