clear all; close all; clc

addpath ./datana-1.1.3/
load ex5_1.mat
%% Plotting batch 1
figure(1);
subplot(1,2,1); plot(1:1:size(cons,1),cons,'-') ;
title('NIR absorbance spectra for batch 1'),xlabel('wavelengths')
subplot(1,2,2); plot(1:1:size(spec,2),spec','- ');
title('concentrations '), xlabel(' time ')

%% Task 2
%PCR of spectra and concerntration
clc
x=1:1000;
t1=1:26;
spec=spec';
figure 
plot(x,spec);
axis tight
xlabel('wavelength'), ylabel('spectra')
title('Wavelength vs spectra')

figure 
plot(t1,cons);
axis tight
xlabel('length'), ylabel('concerntration')
title('length vs concerntration')
%%

 [t,p,r2]=pca(spec');
 plot(t1,r2, '-b*')
 %%
 dim = 8;
[s_fit,b,T,P,q] = pcreg(spec',cons,dim,2);
figure
for k=1:dim;
 subplot(4,2,k)
 cols = [1:4 ] + (k-1)*4;   
 plot(t1',cons,'bo-',t1',s_fit(:,cols),'r-');
 title(['prediction with dimension,' num2str(k)]);
 hold on
end

%%
close all
 clc
load d04_te.dat
process=d04_te;
process=center(process,mean(process));
process = maverage(process,13,0,1);
rmpath  ./datana-1.1.3/


[COEFF, SCORE, LATENT, TSQUARED, EXPLAINED] = pca(process);
r2=cumsum(EXPLAINED);
x1=linspace(0,48,960);
figure
plot(x1,TSQUARED)
xlabel('Time')
ylabel('T^2 PCA')
legend('Process 4')
title('Process 4 after PCA')
