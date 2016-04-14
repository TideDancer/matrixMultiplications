% input A m*n, B n*p
m = 1024;
n = 1024;
p = 1024;

% dense matrix following normal distribution
A = randn(m,n);
B = randn(n,p);
k = log2(n);

% ------------- multiplication routing -------------
T = zeros(n,1);
for i = 1:n
  T(i) = dot(A(:,i), A(:,i)) * dot(B(i,:), B(i,:));
end
[values, iks] = sort(T, 'descend');
J = iks(1:k);

Q = zeros(k,k);
for i = 1:k
  for j = 1:k
    Q(i,j) = dot(A(:,J(i)), A(:,J(j))) * dot(B(J(i),:), B(J(j),:));
  end
end

r = zeros(k,1);
for i = 1:k
  for j = 1:n
    r(i) = r(i) + dot(A(:,J(i)), A(:,j)) * dot(B(J(i),:), B(j,:));
  end
end

w = inv(Q)*r;
AJ = A(:, J);
BJ = B(J, :);
C_approx = AJ*diag(w)*BJ;

% ------------------- compare --------------------------
C = A*B;
error = C - C_approx;
A_norm = norm(A, 'fro');
B_norm = norm(B, 'fro');
C_norm = norm(C, 'fro');
AB_norm = A_norm * B_norm;
error_norm = norm(error, 'fro');

