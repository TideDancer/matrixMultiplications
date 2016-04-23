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
[pdf, cdf] = sample(A,B,'column2norm',1);
ind = getIndexFromPdf(cdf, c, 1);

C = []; R = [];
for i = 1:c
  C = [C  A(:,ind(i))./sqrt(c*pdf(ind(i)))];
  R = [R; B(ind(i),:)./sqrt(c*pdf(ind(i)))];
end

C_approx = C*R; 
  

% ------------------- compare --------------------------
C = A*B;
error = C - C_approx;
A_norm = norm(A, 'fro');
B_norm = norm(B, 'fro');
AB_norm = A_norm * B_norm;
C_norm = norm(C, 'fro');
error_norm = norm(error, 'fro');
