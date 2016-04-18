% input A m*n, B n*p
n = 128;

% dense matrix following normal distribution
A = randn(n,n);
B = randn(n,n);


% ------------- multiplication routing -------------
% build vector v
vec_v = randn(1,n);
V = zeros(n, n^2); % !!!!!!!!!!need sparse representation !!!!!!!!!!!
for i = 1:n
  V(i, n*i-n+1:n*i) = vec_v;
end
% c = zeros(n^2, 1); c is the vector of size n^2 to be solved
u = A*(B*vec_v');
% solve the quation Vc = u


% using steepest descent method as following
% 1. build A = V'V, y =  V'u, then to solve Ac = y, starting by choosing x_i
% r_i = y - A x_i; alpha_i = r_i' r_i / (r_i' A r_i); x_i+1 = x_i + alpha_i * r_i 
F = V'*V;
y = V'*u;
x = randn(n^2, 1); x_last = x.*10;
e = 1e-3; step_limit = n; step = 1;
while norm(x-x_last) > e
  if step > step_limit
    break;
  end
  step = step + 1;
  norm(x-x_last)
  x_last = x;
  r = y - F*x;
  alpha = r'*r / (r'*F*r);
  x = x + alpha * r;
end
% need to be modified !!!!!!!!!!!!!!!!!!!!!!!
% need sparse representation of V, x



% build C from c
C_approx = zeros(n,n);
k = 1;
for i = 1:n
  for j = 1:n
    C_approx(i,j) = x(k); k = k+1;
  end
end


% ------------------- compare --------------------------
C = A*B;
error = C - C_approx;
A_norm = norm(A, 'fro');
B_norm = norm(B, 'fro');
AB_norm = A_norm * B_norm;
C_norm = norm(C, 'fro');
error_norm = norm(error, 'fro');
