% l1 sampling according to the {Achlioptas NearOptimal} paper
% return a sketch B for data matrix A
% input is data matrix A, sampling budget s, failure probability delta
function B = sampleElementL1(A, s, delta);
[m, n] = size(A);

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

B = zeros(m,n);
for k = 1:s
  index = findInInterval(cdf(2:end), rand);
  i = ceil(index/n);
  j = mod(index, n);
  if j == 0
    j = n;
  end
  value = A(i,j)/pdf(index);
  B(i,j) = B(i,j) + value;
end
B = B./s;
return;
