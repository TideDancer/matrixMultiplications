% support several type of J-L based projecting
% parameterList define parameters

function [PA, PB] = project(A, B, type, parameterList);

[ra, ca] = size(A);
[rb, cb] = size(B);
delta = parameterList(1);
epsilon = parameterList(2);
const = parameterList(3);
PA = []; PB = [];

% type == 'kyrillidis2014approximate', return projected matrix PA, PB
% G is gaussian matrix generated following Theorem 1 in the paper
% Prob[||PA*PB-A*B||_2 <= epsilon*A_2*B_2] >= 1-delta
% [delta, epsilon, const] = parameterList, delta and epsilon in (0,1), const is a constant
if strcmp(type, 'kyrillidis2014approximate')
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
  
  % !!!!!!!!!!!!!for approximate nr, set it to small number otherwise too big !!!!!!!!!!!

  t = const*((nra + nrb + log10(log10(1/epsilon)) + log10(1/delta)) / epsilon^2);
  t = round(t);
  G = randn(t,ca).*sqrt(1/t);  % assuming A,B has same size here
  PA = A*G'; PB = G*B;
  return;
end

% classical fast johnson linderstrauss transform
if strcmp(type, 'FJLT')
  PA = fastJLT(A', [delta, epsilon, const]);
  PA = PA';
  PB = fastJLT(B, [delta, epsilon, const]);
end

% according to T.Sarlos 2006 paper algorithm 1, there is tug-of-war projection method achieve same time bound
% return value PA and PB, which can be multiplied as PA*PB to get approximated matrix product
% note that dimension of A or B need to be larger than 16
if strcmp(type, 'tug-of-war')
  k = round(log10(1/delta))
  z = zeros(1,k); SB = zeros(1,k); AS = zeros(1,k);
  for i = 1: k
    S = sign(randn(round(1/epsilon^2), ca));
    SB(i) = S*B;
    AS(i) = A*S'; 
    kk = round(2*(k + log10(k)));
    y = zeros(1,kk);
    for j = 1: kk
      Q = sign(randn(16, cb))
      BQ = B*Q';
      X = A*BQ;
      SBQ = SB*Q';
      X_hat = AS * SBQ;
      y(j) = norm(X-X_hat, 'fro')^2;
    end
    z(i) = median(y);
  end
  [M,i] = min(z);
  PA = AS(i);
  PB = SB(i);
  return; 
end

 
