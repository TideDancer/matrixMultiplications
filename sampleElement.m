% element-wise sampling routing
% include L1 and L2 based routing

function S = sampleElement(A, type, parameterList);
[m, n] = size(A)
	
	
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
  % parameterList is sampling budget s, failure probability delta
  
  s = parameterList(1);
  delta = parameterList(2);
  
  % ------------------ compute row distribution subrouting ------------------
  for i = 1:m
    z(i) = norm(A(i,:), 1); % here ignore constant
  end
  
  alpha = sqrt(log10(m+n)/delta/s);
  beta = log10((m+n)/delta)/(3*s);
  
  % directly use solve() in matlab. but in practice could use binary search because sum(rho) are strictly decreasing order w.r.t. x
  syms x;
  eqn = sum((alpha.*z(1:n)/2/x + sqrt((alpha.*z(1:n)/2/x).^2 + beta.*z(1:n)/x)).^2) == 1;
  solx = double(vpasolve(eqn, x)); % need to use double() to convert sym type to double
  
  rho = zeros(m, 1);
  for i = 1:m
    rho(i) = (alpha*z(i)/2/solx + sqrt((alpha*z(i)/2/solx)^2 + beta*z(i)/solx))^2;
  end
  
  % ------------------ main routing ---------------------
  pdf = zeros(1,m*n);
  cdf = zeros(1,m*n+1);
  for i = 1:m
    for j = 1:n
      pdf((i-1)*n+j) = rho(i)*abs(A(i,j)/z(i));
      cdf((i-1)*n+j+1) = cdf((i-1)*n+j) + pdf((i-1)*n+j); 
    end
  end
  
  S = zeros(m,n);
  for k = 1:s
    index = findInInterval(cdf(2:end), rand);
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
