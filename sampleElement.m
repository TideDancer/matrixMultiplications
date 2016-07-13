% element-wise sampling routing
% return sampled matrix
% include L1 and L2 based routing

function S = sampleElement(A, type, parameterList);
[m, n] = size(A);
	
	
if strcmp(type, 'l2')
  % elementsise sampling based on drines paper, return sampled matrix S
  % sample probability p_ij = min(1, l*A_ij^2 / A_F^2) if |A_ij| > A_F* log^3(2n) / sqrt(2nl)
  % = min(1, sqrt(l)*|A_ij|*log^3(2n) / sqrt(2n) A_F, otherwise
  % l is the parameter, default value set to 1
  % from paper, it's better that l >= 1 and l <= ||A||_F^2 / max A_ij^2
  l = parameterList;
  
  % compute probability
  p = zeros(m, n);
  A_F = norm(A, 'fro');
  A_F2 = A_F^2;
  for i = 1:m
    for j = 1:n
      if abs(A(i,j)) > A_F*log10(2*n)^3
        p(i,j) = min(1, l * A(i,j)^2 / A_F2);
      else
        p(i,j) = min(1, sqrt(l) * abs(A(i,j)) * log10(2*n)^3 / (sqrt(2*n) * A_F) ); 
      end
    end
  end
  
  % do sample
  S = zeros(m,n);
  for i = 1:m
    for j = 1:n
      if rand <= p(i,j)
        S(i,j) = A(i,j)/p(i,j);
      else
        S(i,j) = 0;
      end
    end
  end
  
  return;


elseif strcmp(type, 'l1')
  % l1 sampling according to the {Achlioptas NearOptimal} paper
  % return a sketch S for data matrix A
  % input is data matrix A, 
  % parameterList epsilon, usually set to be 1e-2, failure probability delta, usually set to be 1/n
  
  af = norm(A, 'fro');
  if issparse(A)
    a2 = norm(full(A), 2); % because 2-norm function doesn't work for sparse matrix, convert it to full here
  else
    a2 = norm(A, 2);
  end
  a1 = norm(A, 1);
  sr = af^2/a2^2;
  nd = a1^2/af^2;
  
  epsilon = parameterList(1);
  delta = parameterList(2);
  
  % ------------------ compute row distribution subrouting ------------------
  for i = 1:m
    z(i) = norm(A(i,:), 1); % here ignore constant
  end
  
  % compute s
  nrd = sum(z)/af^2;
  s = nrd * sr / epsilon^2 * log10(n/delta) + sqrt(sr * nd / epsilon^2 * log10(n/delta));
  s = round(s);
  
  alpha = sqrt(log10(m+n)/delta/s);
  beta = log10((m+n)/delta)/(3*s);
  
  % use binary search because sum(rho) are strictly decreasing order w.r.t. mid
  precision = 1e-2;
  % set left to make sure sum > 1, then binary search to find right boundary
  left = alpha * min(z);
  right = m*(left+1); % in case left == 0, put +1 here to avoid right to be set to zero
  while sum((alpha.*z(1:m)/2/right + sqrt((alpha.*z(1:m)/2/right).^2 + beta.*z(1:m)/right)).^2) > 1
    right = 2* right;
  end
  % next with left and right boundary, do binary search
  mid = (left + right) / 2;
  sums = sum((alpha.*z(1:m)/2/mid + sqrt((alpha.*z(1:m)/2/mid).^2 + beta.*z(1:m)/mid)).^2);
  while sums > 1+precision || sums < 1-precision
    if sums > 1+precision
      left = mid;
    else
      right = mid;
    end
    mid = (left + right) /2;
    sums = sum((alpha.*z(1:m)/2/mid + sqrt((alpha.*z(1:m)/2/mid).^2 + beta.*z(1:m)/mid)).^2);
  end

  rho = zeros(m, 1);
  for i = 1:m
    rho(i) = (alpha*z(i)/2/mid + sqrt((alpha*z(i)/2/mid)^2 + beta*z(i)/mid))^2;
  end
  
  % ------------------ main routing ---------------------
  pdf = zeros(1,m*n);
  for i = 1:m
    for j = 1:n
      if z(i) == 0
        pdf((i-1)*n+j) = 0;       
      else
        pdf((i-1)*n+j) = rho(i)*abs(A(i,j)/z(i));
      end
    end
  end
  pdf = pdf./sum(pdf);
  
  S = zeros(m,n);
  % sample s elements
  indices = datasample(1:length(pdf), s, 'Replace', true, 'Weights', pdf);
  for k = 1:s
    index = indices(k);
    i = ceil(index/n);
    j = mod(index, n);
    if j == 0
      j = n;
    end
    value = A(i,j)/pdf(index);
    S(i,j) = S(i,j) + value;
  end
  S = S./s;
  return;

end
