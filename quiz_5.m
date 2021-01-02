%% Task 2
%PCR of spectra and concerntration
clc
load quizz5.mat 
x=1:1000;
t1=1:26;
spec=spec';
figure 
plot(x,spec);
axis tight
xlabel('wavelength'), ylabel('spectra')
title('Wavelength vs spectra')

figure 
plot(t1,c);
axis tight
xlabel('length'), ylabel('concerntration')
title('length vs concerntration')


 [t,p,r2]=pca(spec);
 plot(t1,r2, '-b*')
 
 dim = 25;
[s_fit,b,T,P,q] = pcreg(spec',c,dim,2);

figure()
for k=1:dim;
 cols = [1:4 ] + (k-1)*4;   
 plot(t1',c,'bo-',t1',s_fit(:,cols),'r-');
 title(['prediction with dimension,' num2str(k)]);pause
end
pause
close all
 

%check the fit with dim=1 against the know 'truth'
figure()
plot(t1',s,'bo-',t1',s_fit(:,1:4),'r-');pause

%crossvalidation check yet
[Q2,Q2y,ypred,press,pressy]=crospcrq(spec,sr,1:5,5,1);
Q2y

 