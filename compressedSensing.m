% input A m*n, B n*p
% !!!!!!!!! deal with (C,gamma)-compressibe matrix A and B !!!!!!!!!!!!!!!
n = 1024;

% dense matrix following normal distribution
A = randn(n,n);
B = randn(n,n);

% C-gamma compressible
const = 1;
gamma = 1;% c-gamma compressible, check definition
m = round(log10(n));

% build measurement matrix
rows = int64( (m+1/gamma) * 2^(const * log10(log10(n)^2)) );
M = randn(rows, n);

% reconstruct
% need to invoke l1-magic package (compressed sensing reconstruction package written in matlab) and use l1eq_pd() function
% refer to l1-magic package user mannual
addpath 'l1magic/l1magic/Optimization/';

% following is from l1magic examples except the observation building P
M = orth(M')'; % orthogonalize matrix according to examples in l1magic, a transpose at the end ensure dimension match

% observations
P = (M*A)*B;

% errorbound
epsilon = 1e-3;

% solve the LP
C_approx = [];
for i = 1:n
  % initial guess = min energy
  x0 = M'*P(:,i);
	xp = l1eq_pd(x0, M, [], P(:,i), epsilon);
  C_approx = [C_approx, xp];
end


% ------------------- compare --------------------------
C = A*B;
error = C - C_approx;
A_norm = norm(A, 'fro');
B_norm = norm(B, 'fro');
AB_norm = A_norm * B_norm;
C_norm = norm(C, 'fro');
error_norm = norm(error, 'fro');
