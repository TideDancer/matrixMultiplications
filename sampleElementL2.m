% elementsise sampling based on drines paper, return sampled matrix S
% sample probability p_ij = min(1, l*A_ij^2 / A_F^2) if |A_ij| > A_F* log^3(2n) / sqrt(2nl)
% = min(1, sqrt(l)*|A_ij|*log^3(2n) / sqrt(2n) A_F, otherwise
% from paper, it's better that l >= 1 and l <= ||A||_F^2 / max A_ij^2

function S = sampleElementL2(A, l);
[m, n] = size(A);

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
