clear all
clc
close all
n = 300;
r = 20; % radius (maximal 49)
noise =  randn(n);
[x,y]=meshgrid(-r:r,-r:r);
mask=((x.^2+y.^2)<=r^2);  %(2*r+1)x(2*r+1) bit mask
x = zeros(n,n);
nmin = 50; nmax = 250;
for i=nmin:nmax
    for j=nmin:nmax
        A = noise((i-r):(i+r), (j-r):(j+r));
        x(i,j) = sum(sum(A.*mask));
    end
end
Nr = sum(sum(mask));
x = x(nmin:nmax, nmin:nmax)/Nr;imagesc(x); 
colormap(jet)