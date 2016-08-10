% input A m*n, B n*p
% deal with (C,gamma)-compressibe matrix A and B 

function C_approx = compressedSensing(A, B, parameterList);

[r, n] = size(A);

% C-gamma compressible
const = parameterList(1); % default = 1;
gamma = parameterList(2); % default = 1; c-gamma compressible, check definition
m = round(log10(n));

% build measurement matrix
rows = round( (m+1/gamma) * 2^(const * log10(log10(n))^2) );
if rows > r
  'rows > r'
  return;
end
rows
r
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

return;
