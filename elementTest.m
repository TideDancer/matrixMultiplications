function l = elementTest(A, B, epsilon);

[m,n] = size(A); [n,p] = size(B);

% ------------- validity check -------------------
% from theorem 4 in 2006 drines paper
if 2*n < log10(2*n)^6
  'elementTestFail: 2n < log^6 (2n)'
  return;
end

% ------------- compute l based on epsion -------------
lmax = min(norm(A,'fro')^2/max(max(A.^2)), norm(B,'fro')^2/max(max(B.^2)));
if lmax < 1
  disp('elementTestFail: l < 1');
  return;
end

for l = lmax: -(lmax-1)/100 : 1
  eps = 20*sqrt(n/l) + 100*n/l;
  if eps > epsilon
    break;
  end
end
l = l + (lmax-1)/100;

if l > lmax
  disp(strcat('elementTestFail: epsilon not satisfied in l valid range, best epsilon obtained is',num2str(eps)));
end

return;
