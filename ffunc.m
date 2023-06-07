function fout = ffunc(Q,force)

% function f for Lorenz model
%
% input: Q as a column vector
% output: f(Q) as a column vector
%

N = length(Q);

fout(1)     = -Q(N-1)*Q(N) + Q(N)*Q(2) - Q(1) + force;
fout(2)     = -Q(N)*Q(1) + Q(1)*Q(3) - Q(2) + force;
for i = 3:N-1
    fout(i) = -Q(i-2)*Q(i-1) + Q(i-1)*Q(i+1) - Q(i) + force;
end
fout(N)     = -Q(N-2)*Q(N-1) + Q(N-1)*Q(1) - Q(N) + force;

fout        = fout';

end