% in Drines paper
% A m*n, B n*p
% sample c columns of A and corresponding rows of B
% sampleType corresponds to those in sample.m

function C_approx = basicMatrixMult(A, B, sampleType, parameterList);

[r, n] = size(A);
c = round(log10(n));

% ------------- multiplication routing -------------
[pdf, cdf] = sample(A,B,sampleType,parameterList);
ind = datasample(1:length(pdf), c, 'Replace', false, 'Weights', pdf);
ind = sort(ind); % need to be sorted thus add them in orders to C and R

C = []; R = [];
for i = 1:c
  C = [C  A(:,ind(i))./sqrt(c*pdf(ind(i)))];
  R = [R; B(ind(i),:)./sqrt(c*pdf(ind(i)))];
end

C_approx = C*R; 
  
return;
