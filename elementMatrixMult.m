% in Drines paper
% A m*n, B n*p
% sample c columns of A and corresponding rows of B

m = 1024;
n = 1024;
p = 1024;
c = round(log10(n));

A = randn(m,n);
B = randn(n,p);

% ------------- multiplication routing -------------
% the parameter l is to be tuned, here as default 1, but need to check paper
[pdf, cdf] = sample(A, [], 'elementSquare', [n,1]); % l is a parameter needs to be tuned
S = zeros(m,n);
for i = 1:m
  for j = 1:n
    if rand <= pdf((i-1)*c+j)
      S(i,j) = A(i,j)/pdf((i-1)*c+j);
    else
      S(i,j) = 0;
    end
  end
end

[pdf, cdf] = sample(B, [], 'elementSquare', [n,1]); % l is a parameter needs to be tuned
R = zeros(n,p);
for i = 1:n
  for j = 1:p
    if rand <= pdf((i-1)*p+j)
      R(i,j) = B(i,j)/pdf((i-1)*p+j);
    else
      R(i,j) = 0;
    end
  end
end

C_approx = S*R;


% ------------------- compare --------------------------
C = A*B;
error = C - C_approx;
A_norm = norm(A, 'fro');
B_norm = norm(B, 'fro');
AB_norm = A_norm * B_norm;
C_norm = norm(C, 'fro');
error_norm = norm(error, 'fro');
