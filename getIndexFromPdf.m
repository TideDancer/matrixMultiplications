% given pdf, return [num] sample index using [pdf]
% if isCdf == 1, then this pdf should be a cdf here
% index is in sorted order
% assume pdf is already normalized, meaning sum(pdf) = 1

index = getIndexFromPdf(pdf, num, isCdf);

if isCdf == 1
  cdf = pdf;
else
  cdf = pdf;
  for i = 2:length(cdf)
    cdf(i) = cdf(i-1) + pdf(i);
  end
end

for i = 1:num
  index(i) = findInInterval(cdf, rand);
end
index = sort(index);

return;
