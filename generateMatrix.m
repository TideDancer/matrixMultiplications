% generate A: m*n, B: n*p
m = 1024;
n = 1024;
p = 1024;

% dense matrix following normal distribution
A = randn(m,n);
B = randn(n,p);

% dense matrix non-negative elements, uniform distribution
% A = rand(m,n);
% B = rand(m,n);

% sparse matrix
% in armadillo, use sprandu, sprandn
% sprandu(n_row, n_col, density)
% sp_mat A = sprandu<sp_mat>(100, 200, 0.1);


