function [JHess,B_r,H_r,XB_r,Sr,G_BTQ] = alpha_bounded_BT_Quantities_lorenz(r,N,n_obs,endtime,B,R,M,H,XB)

% if rank is chosen too big, take full model
if r>N
    r = N;
end

step    = endtime/n_obs;

%% compute the shifted system
E       = eig(M);
inds    = abs(E)>1;
if sum(inds)>0
    Emax     = max(abs(E(inds)));
    alpha    = abs(Emax) + 1e-5;
    M = M/alpha;
    H = H/sqrt(alpha);
    B = B/alpha;
end

%% set given prior = Reachability Gramian
Gamma_pr    = B;
L_pr        = chol(B);
A           = M*Gamma_pr+Gamma_pr*M';
if sum(real(eig(A))>0) > 0
    disp('The prior covariance L_pr is not prior-compatible.')
end

G   = H'*(R\H);

%% compute Obs Gramian
Q   = dlyap(M',G)';
% compute a square root factorization of Q
% floating point computation errors induce complex zeros
[L,D] = ldl(Q,'lower');
L_Q = real(L)*real(sqrt(D));

%% balancing with Q 
[V,S,W]     = svd(L_Q'*L_pr); 
S           = S(1:r,1:r);
delQ        = diag(S);
Siginvsqrt  = diag(1./sqrt(delQ));
Sr          = (Siginvsqrt*V(:,1:r)'*L_Q')';
Tr          = L_pr*W(:,1:r)*Siginvsqrt; % balancing transformation
M_r         = Sr'*M*Tr;
H_r         = H*Tr;
B_r         = Sr'*B*Sr;
XB_r        = Sr'*XB;

% Balancing with Q - generate G_BT,H_BT
G_BTQ       = zeros(n_obs*(size(H,1)),N+1);
G_BTQo      = zeros(n_obs*(size(H,1)),N+1);
iter        = expm(M_r*step);
temp        = H_r;
for i = 1:n_obs
    temp                                         = temp*iter;
    temp_proj                                    = temp*Sr';
    G_BTQ((i-1)*(size(H,1))+1:i*(size(H,1)),:)   = temp_proj;
    G_BTQo((i-1)*(size(H,1))+1:i*(size(H,1)),:)  = sqrt(R)\temp_proj;
end
L_prinv           = L_pr\eye(N+1);

% Balancing with Q - compute posterior covariance and mean
R_posinv        = qr([G_BTQo; L_prinv],0);
R_posinv        = triu(R_posinv(1:(N+1),:)); % Pull out upper triangular factor
JHess           = R_posinv' * R_posinv;

end