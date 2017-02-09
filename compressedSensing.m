% input A m*n, B n*p
% deal with (C,gamma)-compressibe matrix A and B 

function C_approx = compressedSensing(A, B, parameterList);

[r, n] = size(A);

% C-gamma compressible
const = parameterList(1); % default = 1;
gamma = parameterList(2); % default = 1; c-gamma compressible, check definition
epsilon = parameterList(3);
sampleSize = parameterList(4);
m = round(log10(n));

% build measurement matrix
if sampleSize == []
  rows = round( (m+1/gamma) * 2^(const * log10(log10(n))^2) );
  if rows > r
    disp('rows > r');
    C_approx = inf;
    return;
  end
else
  rows = sampleSize;
end
M = randn(rows, n);

% reconstruct
% need to invoke l1-magic package (compressed sensing reconstruction package written in matlab) and use l1eq_pd() function
% refer to l1-magic package user mannual
addpath 'l1magic/l1magic/Optimization/';
addpath 'YALL1_v1.4';

% following is from l1magic examples except the observation building P
% M = orth(M')'; % orthogonalize matrix according to examples in l1magic, a transpose at the end ensure dimension match

% observations
P = (M*A)*B;

% solve the LP
C_approx = [];
for i = 1:n
  % ------------ l1 magic ---------
  % % initial guess = min energy
  % x0 = M'*P(:,i);
	% xp = l1eq_pd(x0, M, [], P(:,i), epsilon);
  % -------------------------------
  
  % ------------ yall1 ------------
  opt.tol = epsilon;
  [xp,Out] = yall1(M, P(:,i), opt);
  % -------------------------------

  C_approx = [C_approx, xp];
end

return;
