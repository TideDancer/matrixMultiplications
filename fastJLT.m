% perform fast johnson linderstrauss transform by doing FJLT = P*H*D*x, following 2009 Ailon paper, where set p=2 norm
% project r x c matrix A to k x c matrix PA
% input: matrix A, note that A will be treated as column vector set, thus transform will be operated on each columns
% output: projected matrix PA

function PA = fastJLT(A, parameterList)

[r,c] = size(A);

delta = parameterList(1);
epsilon = parameterList(2);
const = parameterList(3);
k = const * log10(c) / epsilon^2;

% compute Dx, where D is diagonal +1 -1 matrix, which is to randomly assign + - sign to each element of A
D = sign(randn(r,c));
PA = D.*A;

% compute HDx, using matlab fwht function to do Fast Walsh-Hadamard transform, will use 'hadamard' ordering and times r
PA = fwht(PA, r, 'hadamard').*r;

% compute PHDx, P is random generated matrix following certain distribution based on the 2009 Ailon paper
q = (log10(c))^2/r;
if q > 1
  q = 1;
end 
P = rand(k, c);
P = (P < q);
P1 = sqrt(1/q) * randn(k, c);
P = P.*P1;
PA = P * PA;

return;