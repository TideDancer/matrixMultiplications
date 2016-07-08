% in Drines paper
% A m*n, B n*p
% sample c columns of A and corresponding rows of B
% sampleType corresponds to sampleElement function input argument, including 'l1' and 'l2'
% parameterList corresponds to sampleElement function, for 'l2', defult value l is 1

function C_approx = elementMatrixMult(A, B, sampleType, parameterList);

% ------------- multiplication routing -------------
% the parameter l is to be tuned, here as default 1, but need to check paper
l = parameterList(1);
S = sampleElement(A, 'l2', l);
R = sampleElement(B, 'l2', l);

% sparsify representation significantly reduce computing time
S = sparse(S);
R = sparse(R);

C_approx = S*R;

return;
