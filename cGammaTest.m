function [prob, gamma, c] = cGammaTest(AB);

[row, col] = size(AB);
gamma = zeros(col,1);
c = ceil(2*abs(max(max(AB))));

cnt = 0;
for i = 1:col
  v = sort(abs(AB(:,i)), 'descend');

  gmax = log2(c/v(1));
  gmin = log2(c/v(end))/row;
  gamma(i) = min(gmax, gmin);
  for l = 1:row
    if v(l) > c * 2^(-gamma(i) * l)
      cnt = cnt + 1;
    end
  end
end

prob = 1 - cnt / (row*col);
return
