% input a n x d matrix A, output rescaled matrix A to do uniform sampling

function S = sampleUniformIter(A, );
[m, n] = size(A);


% uniform sample n/2 rows from original matrix 
[pdf, cdf] = sample(A, [], 'uniform', 'row')
A_now = A;
[r, c] = size(A_now);
A_half = [];

while r > n*log10(n)
  % uniform sample n/2 rows from  matrix 
	ind = zeros(1,r/2);
  [pdf, cdf] = sample(A_now, [], 'uniform', 'row')
	ind = getIndexFromPdf(cdf, r/2, 1);
	A_half = A_now(ind, :);

  % compute approximate leverage scores of A_now: u_i = a_i * (A'A)^+ * a_i'
	l = zeros(1,r);
	for i = 1:r/2
	  l(i) = A_half(i,:) * pinv(A_half' * A_half) * A_half;
	end 
	l = l./sum(l);
	
	% sample r/2 rows using leverage score as pdf
	ind = getIndexFromPdf(l, r/2, 0);
	A_now = A_now(ind, :);
	[r, c] = size(A_now);
end


