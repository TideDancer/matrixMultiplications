% define dimension when run

delta = 1/dim;  % failure probability
epsilon = 1e-1; % error norm <= epsilon ||A|| ||B||
beta = 1;
const = 1;

density = 'dense';
distribution = 'normal';
matrix = 'high condition';
cond_num = 10^5;

% ------- geometrically distributed singular values ------
% A = gallery('randsvd', dim ,cond_num, 3);
% B = gallery('randsvd', dim ,cond_num, 3);
% disp('gallery(randsvd, dim ,cond_num, 3)');

% % ------- square matrix ------
% A = squareMatrixGen(dim, density, distribution);
% B = squareMatrixGen(dim, density, distribution);

% ------- crazy matrix --------
 A = gallery('sampling', dim);
 B = gallery('chebspec',dim,1);
 disp('gallery(sampling,dim),gallery(chebspec,dim,1)');

% -------------------------- regular -------------------------------------
disp('-----------------regular mult ----------------------');
tic;
C = A*B;
toc;


nnzAB = nnz(C);
A_norm = norm(A, 'fro');
B_norm = norm(B, 'fro');
AB_norm = A_norm * B_norm;
C_norm = norm(C, 'fro');

