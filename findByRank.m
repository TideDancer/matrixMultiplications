% find values of rankStart to rankEnd in a 2-D matrix, and return them as a list
% value is the array of those ranked values
% position is the matrix of corresponding locations (i,j) in each row

function [value, position] = findByRank(M, rankStart, rankEnd);

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
  value = [value; maxv];
  position = [position; x y];
  M(x,y) = minv;
end

