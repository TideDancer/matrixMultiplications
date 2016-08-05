% input A m*n, B n*p


function C_approx = nystrom(A, B);

[m,n] = size(A);
k = round(log10(n));

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

return;
