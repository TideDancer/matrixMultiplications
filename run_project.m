% random matrix test 
delta = 1e-2;  % failure probability
epsilon = 1e-2; % error norm <= epsilon ||A|| ||B||
beta = 1;
const = 1;

matrix = 'high condition';
cond_num = 10^5;

r = 2^15;
c = 2^15;

% % build coherent matrix
% Z = zeros(dim); I = eye(dim); O = ones(dim).*1e-8; 
% R = rand(dim).*1e-8; alphaB = randn(dim)*1e8;
% A = [A(1:dim/2, :); Z(1:dim/2, :)] + O;
% B = [B(:, 1:dim/2)  Z(:, 1:dim/2)] + O;

%%%% dexter %%%%%
% load('dexter.mat');
% A = dexter1; B = dexter2';
% min(leverage(A))

% %%%% heart1 %%%%%
% load('goodwin.mat');
% A = full(Problem.A)./1000;
% B = A(randperm(end),randperm(end));

% ------- geometrically distributed singular values ------
% A = gallery('randsvd', dim ,cond_num, 3);
% B = gallery('randsvd', dim ,cond_num, 3);
% disp('gallery(randsvd, dim ,cond_num, 3)');

%% ------- crazy matrix --------
% A = gallery('sampling', dim);
% B = gallery('chebspec',dim,1);
% disp('gallery(sampling,dim),gallery(chebspec,dim,1)');

% start computing
for sampleDim = 5:14
sampleSize = 2^sampleDim;
disp(sampleSize);
disp('--------------------------------------------');
disp('--------------------------------------------');
for k = 1:5
  % ------- randn matrix ------
  clear A;
  clear B;
  A = randn(r,c);
  B = randn(r,c);

  disp('generating done');
  disp('direct mult');
  tic;
  C = A*B;
  toc;
  disp('-----------------------');
  A_norm = norm(A, 'fro');
  B_norm = norm(B, 'fro');
  AB_norm = A_norm * B_norm
  C_norm = norm(C, 'fro')

  for i = 1:10
    disp('--------- following methods -----------');
    tic;
    C_approx = randomProjMult(A, B, 'FJLT', [delta, epsilon, const, sampleSize]);
    toc;
    error = C - C_approx;
    error_norm = norm(error, 'fro');
    disp(error1_norm/AB_norm);

    tic;
    C_approx = randomProjMult(A, B, 'tug-of-war', [delta, epsilon, const, sampleSize]);
    toc;
    error = C - C_approx;
    error_norm = norm(error, 'fro');
    disp(error1_norm/AB_norm);
    
    tic;
    C_approx = randomProjMult(A, B, 'clarkson2009numerical', [delta, epsilon, const, sampleSize]);
    toc;
    error = C - C_approx;
    error_norm = norm(error, 'fro');
    disp(error1_norm/AB_norm);

    tic;
    C_approx = randomProjMult(A, B, 'kyrillidis2014approximate', [delta, epsilon, const, sampleSize]);
    toc;
    error = C - C_approx;
    error_norm = norm(error, 'fro');
    disp(error1_norm/AB_norm);
  end
end
end
