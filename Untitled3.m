n = 2^8;
t1 = [0:1/n:1-1/n];
t2 = t1;
for i=1:n   % first row of cov. matrix, arranged in a matrix
    for j=1:n
        G(i,j)=exp(-8*sqrt(min(abs(t1(1)-t1(i)), ...
            1-abs(t1(1)-t1(i)))^2 + min(abs(t2(1)-t2(j)), ...
            1-abs(t2(1)-t2(j)))^2));
    end
end
Gamma = fft2(G);   % the eigenvalue matrix n*fft2(G/n)
Z = randn(n,n) + sqrt(-1)*randn(n,n);
X = real(fft2(sqrt(Gamma).*Z/n));
imagesc(X); 
colormap(jet)