%%
clear all
close all
clc

addpath('C:\Users\UTHMAN\Documents\Advanced_DA&ML/drtoolbox')

load candy_matlab.mat

no_dims = intrinsic_dim(candy)
normal=normalize(candy,'norm');

no_dims1 = intrinsic_dim(normal)
[coeff, score]=pca(normal);

rng('default')
[W,H, D] = nnmf(normal,2);

%D is the root mean square error
%recX = reconstruct_data(H, coeff)

%error=sqrt(mse(recX,normal))

