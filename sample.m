% support several types of sampling
% parameterList define parameters
% if no parameters, set this argument to be [] then default parameters will be set
% pdf return probability density of each elements
% cdf return cummulative density

function [pdf, cdf] = sample(A, B, types, parameterList);

[ra, ca] = size(A);
[rb, cb] = size(B);

% if type == 'column2norm', return p(i) = beta * A_k_column_2norm * B_k_row_2norm / sum_j(A_j_column_2norm * B_j_row_2norm)
% beta = parameterList
% don't put parameterList to be more than one value, don't put parameterList <= 0
if types == 'column2norm'
  beta = parameterList;
  p = zeros(1, ca);
  p_sum = 0;
  if ca ~= rb
    return;
  end 
  for i = 1:ca
    p(i) = norm(A(:,i)) * norm(B(i,:)) * beta;
    p_sum = p_sum + p(i);
  end
  for i = 1:ca
    p(i) = p(i)/p_sum;
  end
end


% -------- return pdf and cdf ---------
if length(p) > 0
  pdf = p;
  cdf = pdf;
  for i = 2:length(pdf)
    cdf(i) = cdf(i-1) + pdf(i);
  end
end
return
