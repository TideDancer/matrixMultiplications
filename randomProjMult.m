% matrix multiplication using random projection
% input A m*n, B n*p
% there are several project methods

function C_approx = randomProjMult(A, B, type, parameterList);

% do projection based on the paper kyrillidis2014approximate
% [PA, PB] = project(A,B, 'kyrillidis2014approximate',[1e-3, 1e-3, 1e-6]);
% this corresponds to kyrillidis2014 paper, also other project method like 'fjlt', 'tug of war' can be used

[PA, PB] = project(A, B, type, parameterList);

C_approx = PA*PB;

return;
