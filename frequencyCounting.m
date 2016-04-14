% very problematic !!!!!!!!!!!!!!!!!!!!!!!!!!!!

% input A n*n, B n*n, A,B are non-negative matrix
n = 256;
A = rand(n,n);
B = rand(n,n);
b = log2(n);

% ------------- compute summary -------------
S = zeros(b,3);
for i = 1:n
  R = A(:,i)*B(i,:);
  tmp = findByRank(R, b+1, b+1);
  w = tmp(1,1);
  
  L = findByRank(R, 1, b);
  L(:,1) = L(:,1) - w;
  
  % for each entry e occuring in both S and L, double e's weight in S and remove e from L
  % get S union L
  SIL = intersect(S, L, 'rows');
  S_left = setdiff(S, SIL, 'rows');
  L_left = setdiff(L, SIL, 'rows');
  SIL(:,1) = 2.*SIL(:,1);
  SUL = flipud(union(union(S_left, L_left, 'rows'), SIL, 'rows', 'sorted'));
  
  [r,c] = size(SUL);
  if r > b
  	w = SUL(b+1, 1);
  else
    w = 0;
  end
  S = SUL(1:b,:);
  S(:,1) = S(:,1) - w;
end
 
% --------------- get approx value for each entry -----------
C_approx = zeros(n,n);
for i = 1:n
  for j = 1:n
    [lia, loc] = ismember([i j], S(:, 2:3), 'rows');
    if lia == 1
      C_approx(i,j) = S(loc, 1);
    end
  end
end



% ------------------- compare --------------------------
C = A*B;
error = C - C_approx;
A_norm = norm(A, 'fro');
B_norm = norm(B, 'fro');
C_norm = norm(C, 'fro');
AB_norm = A_norm * B_norm;

error_norm = norm(error, 'fro');
