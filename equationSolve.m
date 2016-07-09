% input A m*n, B n*p

function C_approx = equationSolve(A, B);

[r, n] = size(A);

% ------------- multiplication routing -------------
% build vector v
vec_v = randn(1,n);

% V' is not necessary, so commented here
% V = sparse(n, n^2); 
% for i = 1:n
%   V(i, n*i-n+1:n*i) = vec_v;
% end

% c = zeros(n^2, 1); c is the vector of size n^2 to be solved
u = A*(B*vec_v');
% solve the quation Vc = u


% using steepest descent method as following
% 1. build A = V'V, y =  V'u, then to solve Ac = y, starting by choosing x_i
% r_i = y - A x_i; alpha_i = r_i' r_i / (r_i' A r_i); x_i+1 = x_i + alpha_i * r_i 

% compute F = V'*V, F is block diagonal, blocks are the same, in total n block members;
F_block = vec_v'*vec_v; F_num = n;

% compute y = V'*u;
y = [];
for i = 1:n
  y = [y; vec_v'.*u(i)];
end
  
x = randn(n^2, 1); x_last = x.*10;
e = 1e-3; step_limit = n; step = 1;
while norm(x-x_last) > e
  step
  if step > step_limit
    break;
  end
  step = step + 1;
  x_last = x;
  %r = y - F*x;
  fx = [];
  for i = 1:n
    fx = [fx; F_block * x((i-1)*n+1:i*n)];
  end
  r = y - fx;

  %compute F*r
  fr = [];
  for i = 1:n
    fr = [fr; F_block * r((i-1)*n+1:i*n)];
  end

  alpha = r'*r / (r'*fr);
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

return;
