% support several type of J-L based projecting
% parameterList define parameters

function [PA, PB] = project(A, B, type, parameterList);

[ra, ca] = size(A);
[rb, cb] = size(B);
PA = []; PB = [];

% type == 'kyrillidis2014approximate', return projected matrix PA, PB
% G is gaussian matrix generated following Theorem 1 in the paper
% Prob[||PA*PB-A*B||_2 <= epsilon*A_2*B_2] >= 1-delta
% [delta, epsilon, const] = parameterList, delta and epsilon in (0,1), const is a constant
if strcmp(type, 'kyrillidis2014approximate')
  [delta, epsilon, const] = parameterList;
  if ca ~= rb
    return;
  end
  [ua, sa, va] = svd(A,'econ');
  [ub, sb, vb] = svd(B,'econ');
  % nuclear rank defined as nuclear_norm/2norm
  % nuclear norm is sum of singular values, 2norm is largest singular value
  % !!!!!!!!!!!!! NOTE THAT THIS STEP OF SVD TAKES TOO MUCH TIME AS IT'S AS HARD AS ORIGINAL PROBLEM !!!!!!!!!!!!!!!!!
  nra = sum(max(sa))/max(max(sa));
  nrb = sum(max(sb))/max(max(sb));

  t = const*((nra + nrb + log2(log2(1/epsilon)) + log2(1/delta)) / epsilon^2);
  G = randn(t,d).*sqrt(1/t);
  PA = A*G'; PB = G*B;
  return;
end


