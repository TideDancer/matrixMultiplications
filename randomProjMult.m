% matrix multiplication using random projection
% input A m*n, B n*p
% there are several project methods

function C_approx = randomProjMult(A, B, type, parameterList);

delta = parameterList(1);
const = parameterList(3);
sampleSize = parameterList(4);

if sampleSize == []
  [PA, PB] = project(A, B, type, parameterList);
else

  [ra, ca] = size(A);
  [rb, cb] = size(B);
  PA = []; PB = [];
  
  % type == 'kyrillidis2014approximate', return projected matrix PA, PB
  % G is gaussian matrix generated following Theorem 1 in the paper
  % Prob[||PA*PB-A*B||_2 <= epsilon*A_2*B_2] >= 1-delta
  % [delta, epsilon, const] = parameterList, delta and epsilon in (0,1), const is a constant
  if strcmp(type, 'kyrillidis2014approximate')
    t = sampleSize;
    G = randn(t,ca).*sqrt(1/t);  % assuming A,B has same size here
    PA = A*G'; PB = G*B;
  end
  
  % classical fast johnson linderstrauss transform
  if strcmp(type, 'FJLT')
    PA = fastJLT(A', parameterList);
    PA = PA';
    PB = fastJLT(B, parameterList);
  end
  
  % in 2009 woodruff paper:
  if strcmp(type, 'clarkson2009numerical')
    m = sampleSize;
    S = sign(randn(ca, m));
    PA = A*S./sqrt(m);
    PB = S'*B./sqrt(m);
  end
  
  % according to T.Sarlos 2006 paper algorithm 1, there is tug-of-war projection method achieve same time bound
  % return value PA and PB, which can be multiplied as PA*PB to get approximated matrix product
  % note that dimension of A or B need to be larger than 16
  if strcmp(type, 'tug-of-war')
    k = round(log10(1/delta));
    kk = round(2*(k + log10(k)));
    rs = sampleSize;
    z = zeros(1,k); SB = zeros(rs, cb, k); AS = zeros(ra, rs, k);
    for i = 1: k
      S = sign(randn(rs, ca));
      SB(:,:,i) = S*B;
      AS(:,:,i) = A*S'; 
      y = zeros(1,kk);
      for j = 1: kk
        Q = sign(randn(16, cb));
        BQ = B*Q';
        X = A*BQ;
        SBQ = SB(:,:,i)*Q';
        X_hat = AS(:,:,i) * SBQ;
        y(j) = norm(X-X_hat, 'fro')^2;
      end
      z(i) = median(y);
    end
    [M,i] = min(z);
    PA = AS(:,:,i);
    PB = SB(:,:,i);
  end

end

C_approx = PA*PB;

return;
