% select algorithm from DRINEAS06 paper
% not sure how to use it

D = 0;
for i = 1:n
  D = D + a(i);
  p = a(i)/D;
  i_star = i;
  a_star(i) = a(i);
end
