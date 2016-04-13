% input A m*n, B n*p
n = 1024;

% dense matrix following normal distribution
A = randn(n,n);
B = randn(n,n);

% C-gamma compressible
const = 10;
gamma = % c-gamma compressible, check definition

% build sensing matrix
rows = (m+1/gamma) * 2^(const * log2(log2(n)^2));
M = randn(rows, n);
P = (M*A)*B;

% reconstruct


