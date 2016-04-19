% support several type of sampling
% parameterList define parameters
% if no parameters, set this argument to be [] then default parameters will be set
% pdf return probability density of each elements
% cdf return cummulative density

function [pdf, cdf] = sample(A, B, type, parameterList);

[ra, ca] = size(A);
[rb, cb] = size(B);

% type == 'column2norm', return p(i) = beta * A_k_column_2norm * B_k_row_2norm / sum_j(A_j_column_2norm * B_j_row_2norm)
% beta = parameterList
% don't put parameterList to be more than one value, don't put parameterList <= 0
if strcmp(type, 'column2norm')
  beta = parameterList;
  p = zeros(1, ca);
  p_sum = 0;
  if ca ~= rb
    return;
  end 
  for i = 1:ca
    p(i) = norm(A(:,i)) * norm(B(i,:)) * beta;
    p_sum = p_sum + p(i)/beta;
  end
  for i = 1:ca
    p(i) = p(i)/p_sum;
  end
end


% type == 'column2normFro', return p(i) = beta * A_k_column_2norm * B_k_row_2norm / (A_fro*B_fro)
% beta = parameterList
% don't put parameterList to be more than one value, don't put parameterList <= 0
if strcmp(type, 'column2normFro')
  beta = parameterList;
  p = zeros(1, ca);
  p_sum = norm(A, 'fro') * norm(B, 'fro');
  if ca ~= rb
    return;
  end 
  for i = 1:ca
    p(i) = norm(A(:,i)) * norm(B(i,:)) * beta / p_sum;
  end
end


% type == 'elementSquare', return elementwise sampling probability p(i,j), but in a list format
% finally in the list p(i,j) = p((i-1)*c+j) 
% based on drines paper
% only accept one matrix A and compute p corresponding to A, thus B can be anything
% parameter list is [n, l], where n is inner dimension
if strcmp(type, 'elementSquare')
  if length(parameterList) ~= 2
    disp('put [n, l] as parameterList');
    return;
  end
  n = parameterList(1); l = parameterList(2);
  p = zeros(1,ra*ca);
  A_F = norm(A, 'fro');
  A_F2 = A_F^2;
  for i = 1:ra
    for j = 1:ca
      if abs(A(i,j)) > A_F*log10(2*n)^3
        p((i-1)*ca+j) = min(1, l * A(i,j)^2 / A_F2);
      else
        p((i-1)*ca+j) = min(1, sqrt(l) * abs(A(i,j)) * log10(2*n)^3 / (sqrt(2*n) * A_F) ); 
      end
    end
  end
end


% type = 'column2nomrSquareFro', return p(i) = beta * A_k_column_2norm^2 / A_frob^2
% beta = parameterList
% don't put parameterList to be more than one value, don't put parameterList <= 0
% only matrix A will be used
if strcmp(type, 'column2normSquareFro')
  beta = parameterList;
  p = zeros(1, ca);
  p_sum = norm(A, 'fro');
  for i = 1:ca
    p(i) = norm(A(:,i))^2 * beta / p_sum;
  end
end


% type = 'column2nomrSquare2norm', return p(i) = beta * A_k_column_2norm^2 / sum_j(A_j_column_2norm^2)
% beta = parameterList
% don't put parameterList to be more than one value, don't put parameterList <= 0
% only matrix A will be used
if strcmp(type, 'column2normSquare2norm')
  beta = parameterList;
  p = zeros(1, ca);
  p_sum = 0;
  for i = 1:ca
    p(i) = norm(A(:,i))^2 * beta;
    p_sum = p_sum + p(i)/beta;
  end
  for i = 1:ca
    p(i) = p(i)/p_sum;
  end
end


% type = 'uniform', will return uniform sampling p(i) = 1/n
% isRow = parameterList
% unless isRow == 'row', will do column sample
% only matrix A will be used
if strcmp(type, 'uniform')
  isRow = parameterList;
  if strcmp(isRow, 'row')
    p = zeros(1, ra);
    n = ra;
  else
    p = zeros(1, ca);
    n = ca;
  end
  p_sum = 0;
  for i = 1:n
    p(i) = 1/n;
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
