% input A n*n, B n*n, A,B are non-negative matrix
n = 1024;
A = rand(n,n);
B = rand(n,n);
b = log2(n);

% ------------- compute summary -------------
S = [];
for i = 1:n
  R = A(:,i)*B(i,:);
  [w tmp] = findMax(R, b+1, b+1);
  
  [L_weight, L_position] = findMax(R, 1, b);
  L_weight = L_weight - w;
  
  for j = 1:len(S_weight)
    if ismember(S_weight(j), L_weight)
      S_weight(j) = S_weight(j) + S_weight(j);
    end
  end
  
  
    



% ------------------- compare --------------------------
C = A*B;
error = C - C_approx;
A_norm = norm(A, 'fro');
B_nomr = norm(B, 'fro');
C_norm = norm(C, 'fro');
error_norm = norm(error, 'fro');
