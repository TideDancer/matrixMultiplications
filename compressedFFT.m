% product AB need to has at most b/8 nonzero entries
% d >= 6 log n
% then correctness probability is 1-o(1)

function C_approx = compressedFFT(A, B, nnzAB);

[r, n] = size(A);
b = nnzAB * 8;
if b > n*n
  disp('compressedFFT test fail: b > n*n');
  C_approx = inf;
  return;
end

d = 6*log10(n);

% --------------------- compressed product --------------
s1_t = []; s2_t = []; h1_t = []; h2_t = []; p_t = [];
for t = 1:d
    s1 = sign(randn(1,n));
    s2 = sign(randn(1,n));
    h1 = ceil(rand(1,n).*b);
    h2 = ceil(rand(1,n).*b);
    p = zeros(1,b);
    
    for k = 1:n
        pa = zeros(1,b); pb = zeros(1,b);
        for i = 1:n
            pa(h1(i)) = pa(h1(i))+s1(i)*A(i,k);
        end
        for j = 1:n
            pb(h2(j)) = pb(h2(j))+s2(j)*B(k,j); % original paper: pa(h2(j)) = pb(h2(j)) + s2(j)*B_kj
        end
        pa = fft(pa); pb = fft(pb);
        p = p + pa.*pb;
    end
    
    s1_t = [s1_t; s1];
    s2_t = [s2_t; s2];
    h1_t = [h1_t; h1];
    h2_t = [h2_t; h2];
    p_t = [p_t; p];
end

for t = 1:d
    p_t(t,:) = ifft(p_t(t,:));
end

% ------------------- decompress (i,j) ------------------
C_approx = zeros(n,n);
for i = 1:n
    for j = 1:n
        for t = 1:d
            X(t) = s1_t(t,i)*s2_t(t,j)* p_t(t, 1+mod((h1_t(t,i)+h2_t(t,j)), b));
        end
        C_approx(i,j) = median(X);
    end
end

return;
