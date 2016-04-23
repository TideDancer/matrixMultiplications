% find values of rankStart to rankEnd in a 2-D matrix, and return them as a table
% the first column is the array of those ranked values
% the next two columns are the corresponding locations (i,j) 

function table = findByRank(M, rankStart, rankEnd);

table = [];
value = [];
position = [];

if rankStart <= 0 || rankStart > rankEnd
  return;
end

minv = min(M(:)) - 1;

for i = 1:rankStart-1
  [maxv idx] = max(M(:));
  [x y] = ind2sub(size(M), idx);  
  M(x,y) = minv;
end

for i = rankStart:rankEnd
  [maxv idx] = max(M(:));
  [x y] = ind2sub(size(M), idx);  
  table = [table; maxv x y];
  M(x,y) = minv;
end

