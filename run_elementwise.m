% ############# we can use ||C - C_approx || as a measurement of elementwise error #############################

% random matrix test 
A_norm = [];
B_norm = [];
C_norm = [];
error1_norm = [];
error2_norm = [];

for i = 10:15
  dim = 2^i;
  
  A = squareMatrixGen(dim, 'dense', 'normal');
  B = squareMatrixGen(dim, 'dense', 'normal');
  
  l = min(norm(A, 'fro')^2/max(max(A))^2, norm(B, 'fro')^2/max(max(B))^2);
  epsilon = 1e-2;
  
  C_approx_l2 = elementMatrixMult(A, B, 'l2', l); % default failure probability 1/n
  C_approx_l1 = elementMatrixMult(A, B, 'l1', [epsilon, 1/dim]); % failure probablity 1/n
  
  % ------------------- compare --------------------------
  C = A*B;
  error_2 = C - C_approx_l2;
  error_1 = C - C_approx_l1;
  A_norm = [A_norm; norm(A, 'fro')];
  B_norm = [B_norm; norm(B, 'fro')];
  C_norm = [C_norm; norm(C, 'fro')];
  error2_norm = [error2_norm; norm(error_2, 'fro')];
  error1_norm = [error1_norm; norm(error_1, 'fro')];
end

i = 10:15;
save('result_element_dense_normal.mat', i, A_norm, B_norm, C_norm, error1_norm, error2_norm);

