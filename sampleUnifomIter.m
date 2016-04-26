% input a n x d matrix A, use uniform sampling to estimate leverage score then sample rows based on leverage scores
% apply to matrix that is thin, i.e. n >> d
% return spectral approximation S of A
% refer to 'Uniform Sampling for Matrix Approximation, M.B. Cohen'

function S = sampleUniformIter(A );
[m, n] = size(A);


% uniform sample n/2 rows from original matrix 
A_now = A;
[r, c] = size(A_now);

while r > n*log10(n)
  % uniform sample n/2 rows from  matrix 
  r1 = ceil(r/2);
	ind = zeros(1,r1);
  [pdf, cdf] = sample(A_now, [], 'uniform', 'row');
	ind = findIndexFromPdf(cdf, r1, 1);
	A_half = A_now(ind, :);

  % compute approximate leverage scores of A_now: u_i = a_i * (A'A)^+ * a_i'
	l = zeros(1,r1);
	for i = 1:r1
	  l(i) = A_half(i,:) * pinv(A_half' * A_half) * A_half(i,:)';
	end 
	l = l./sum(l);
	
	% sample r/2 rows using leverage score as pdf
	ind = findIndexFromPdf(l, r1, 0);
	A_now = A_now(ind, :);
	[r, c] = size(A_now);
end

S = A_now;
return

