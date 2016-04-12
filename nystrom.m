% input A m*n, B n*p
k = log2(n);

% ------------- multiplication routing -------------
T = zeros(n);
for i = 1:n
  T(i) = dot(A(i,:), A(i,:)) * dot(B(:,i), B(:,i));
end
[values, iks] = sort(T, 'descend');
J = iks(1:k);



% ------------------- compare --------------------------
C = A*B;
error = C - C_approx;
A_norm = norm(A, 'fro');
B_nomr = norm(B, 'fro');
C_norm = norm(C, 'fro');
error_norm = norm(error, 'fro');
