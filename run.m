% run_time = 50;
% er = zeros(run_time, 8);
% 
% for i = 1:run_time
%   basicMatrixMult;
%   er(i,1) = error_norm/AB_norm;
% 
%   elementMatrixMult;
%   er(i,2) = error_norm/AB_norm;
% 
%   randomProjMult;
%   er(i,3) = error_norm/AB_norm;
% 
%   nystrom;
%   er(i,4) = error_norm/AB_norm;
%   
%   compressedSensing;
%   er(i,5) = error_norm/AB_norm;
%   
%   equationSolve;
%   er(i,6) = error_norm/AB_norm;
% 
%   frequencyCounting;
%   er(i,7) = error_norm/AB_norm;
%   
%   compressedFFT;
%   er(i,8) = error_norm/AB_norm;
% end
% 
% for j = 1:8
%   er(i+1,j) = sum(er(1:i,j));
% end


% ############# we can use ||C - C_approx || as a measurement of elementwise error #############################

% random matrix test 
dim = 2^10;

A = squareMatrixGen(dim, 'dense', 'normal');
B = squareMatrixGen(dim, 'dense', 'normal');

C_approx = basicMatrixMult(A, B, 'column2norm', 1); % coloum 2-norm based sampling

C_approx = elementMatrixMult(A, B, 'l2', 1); % l-2 based element-wise sampling

C_approx = randomProjMult(A, B, 'kyrillidis2014approximate', [1e-3, 1e-3, 1e-6]); % kyrillidis2014 paper

C_approx = nystrom(A, B);




% ------------------- compare --------------------------
C = A*B;
error = C - C_approx;
A_norm = norm(A, 'fro');
B_norm = norm(B, 'fro');
AB_norm = A_norm * B_norm;
C_norm = norm(C, 'fro');
error_norm = norm(error, 'fro');

