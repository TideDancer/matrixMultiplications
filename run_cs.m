settings;

% ------------------------- cs based ---------------------------
disp('--------------compressed sensing---------------');
b = log(nnzAB)/log(dim) < 0.29462; % corollary 5 in the paper of cs-based mult
[prob, gamma, c] = cGammaTest(C);  % c-gamma test to compute gamma for a specific C
tic;
C_approx = compressedSensing(A, B, [1, min(gamma), epsilon]);
toc;
disp(strcat('errorRatio = ', num2str(norm(C-C_approx, 'fro')/AB_norm)))


