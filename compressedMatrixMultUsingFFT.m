n = 2048;
b = log2(n);
d = 20;
A = randn(n,n); B = randn(n,n);

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
            pb(h2(j)) = pb(h2(j))+s2(j)*B(j,k);
        end
        pa = fft(pa); pb = fft(pb);
        for z = 1:b
            p(z) = p(z)+pa(z)*pb(z);
        end
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
        C_approx(i,j) = mean(X);
    end
end

% ------------------- compare --------------------------
C = A*B;
error = C - C_approx;
A_norm = norm(A, 'fro');
B_norm = norm(B, 'fro');
AB_norm = A_norm * B_norm;
C_norm = norm(C, 'fro');
error_norm = norm(error, 'fro');
