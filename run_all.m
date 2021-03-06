settings;

% ------------------------ column l2 sample -----------------------------------
disp('-------------column 2 norm sampling-----------');
tic;
C_approx = basicMatrixMult(A, B, 'column2norm', [delta, epsilon, beta]);
toc; 
disp(strcat('errorRatio = ', num2str(norm(C-C_approx, 'fro')/AB_norm)))


% ------------------------ elementwise l2 ------------------------------------
disp('-------------elementwise-------------');
l = elementTest(A,B,epsilon); % default delta is 1/dim
tic;
C_approx = elementMatrixMult(A, B, 'l2', l); % l-2 based element-wise sampling
toc;
disp(strcat('errorRatio = ', num2str(norm(C-C_approx, 'fro')/AB_norm)))


% ------------------------ fjlt ----------------------------------------
disp('--------------fjlt -----------');
tic;
C_approx = randomProjMult(A, B, 'FJLT', [delta, epsilon, const]);
toc;
disp(strcat('errorRatio = ', num2str(norm(C-C_approx, 'fro')/AB_norm)))


% ------------------------ tug-of-war ----------------------------------
disp('--------------tug-of-war-------------');
tic;
C_approx = randomProjMult(A, B, 'tug-of-war', [delta, epsilon, const]);
toc;
disp(strcat('errorRatio = ', num2str(norm(C-C_approx, 'fro')/AB_norm)))


% -------------------------- clarkson2009numerical -----------------------------
disp('------------clarkson2009numerical---------------');
tic;
C_approx = randomProjMult(A, B, 'clarkson2009numerical', [delta, epsilon, const]); 
toc;
disp(strcat('errorRatio = ', num2str(norm(C-C_approx, 'fro')/AB_norm)))


% ------------------------- kyrillidis2014approximate ---------------------------
disp('--------------kyrillidis2014approximate--------------');
tic;
C_approx = randomProjMult(A, B, 'kyrillidis2014approximate', [delta, epsilon, const]);
toc;
disp(strcat('errorRatio = ', num2str(norm(C-C_approx, 'fro')/AB_norm)))


% --------------------------- compressedFFT ---------------------------
disp('--------------compressedFFT---------------');
tic;
C_approx = compressedFFT(A, B, nnzAB);
toc;
disp(strcat('errorRatio = ', num2str(norm(C-C_approx, 'fro')/AB_norm)))


% --------------------------- equationSolve ---------------------------
disp('-------------equationSolve---------------');
tic;
C_approx = equationSolve(A, B);
toc;
disp(strcat('errorRatio = ', num2str(norm(C-C_approx, 'fro')/AB_norm)))


% -------------------------- nystrom --------------------------------------------
disp('--------------nystrom---------------');
tic;
C_approx = nystrom(A, B, round(log10(dim)));
toc;
disp(strcat('errorRatio = ', num2str(norm(C-C_approx, 'fro')/AB_norm)))


% ------------------------- cs based ---------------------------
disp('--------------compressed sensing---------------');
b = log(nnzAB)/log(dim) < 0.29462; % corollary 5 in the paper of cs-based mult
[prob, gamma, c] = cGammaTest(C);  % c-gamma test to compute gamma for a specific C
tic;
C_approx = compressedSensing(A, B, [1, min(gamma), epsilon]);
toc;
disp(strcat('errorRatio = ', num2str(norm(C-C_approx, 'fro')/AB_norm)))



% % --------------------------- frequencyCounting -----------------------
% tic; 
% C_approx = frequencyCounting(A, B); % A, B need to be non-negative matrix !!!!!!!!!!
% toc;
% disp(strcat('errorRatio = ', num2str(norm(C-C_approx, 'fro')/AB_norm)))
% disp('------------------------------');



