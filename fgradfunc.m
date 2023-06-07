function fgradout = fgradfunc(Q)

% function df/dq for Lorenz model
%
% input: Q (as a row vector) 
% output: df/dq (Q) 
%

N           = length(Q);

A           = -sparse(diag(ones(N,1)));

A(1,2)      = Q(N);
A(1,N-1)    = -Q(N);
A(1,N)      = Q(2) - Q(N-1);
A(2,1)      = -Q(N) + Q(3);
A(2,3)      = Q(1);
A(2,N)      = -Q(1);

for i = 3:N-1
    A(i,i+1) = Q(i-1);
    A(i,i-1) = -Q(i-2) + Q(i+1);
    A(i,i-2) = -Q(i-1);
end

A(N,1)      = Q(N-1);
A(N,N-2)    = -Q(N-1);
A(N,N-1)    = -Q(N-2) + Q(1);
    
fgradout    = A;       

end
