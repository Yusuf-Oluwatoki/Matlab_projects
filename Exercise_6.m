%% Task 1
% PLS regression and PCA regression
clc
load data_5_1.mat

y=c;
%for i=1:29
    [xl,yl,xs,ys,beta,pctvar,mse] = plsregress(X,y,3);
    %plot(1:10,cumsum(100*pctvar(2,:)),'-bo');
    figure
    plot(1:3,cumsum(100*pctvar(2,:)),'-bo');
    yFitted = [ones(size(X,1),1) X]*beta;
    figure
    plot(y,yFitted,'o');
%end

[PCALoadings,PCAScores,PCAVar] = pca(X,'Economy',false);
betaPCR = regress(y-mean(y), PCAScores(:,1:2));
plot(1:3,cumsum(100*pctvar(2,:)),'-bo')

betaPCR = PCALoadings(:,1:2)*betaPCR;
betaPCR = [mean(y) - mean(X)*betaPCR; betaPCR];
yfitPCR = [ones(size(X,1),1) X]*betaPCR;
plot(y,yfitPCR,'o');

figure
plot(y,yFitted,'bo',y,yfitPCR,'r^');
xlabel('Observed Response');
ylabel('Fitted Response');
legend({'PLSR with 3 Components' 'PCR with 3 Components'},  ...
	'location','NW');
 
%% Task 2
 %PCA of spectra and concerntration
clc
load nirdemo.mat 
x=1:251;

figure 
subplot(2,2,1)
plot(x,x1);
axis tight
xlabel('wavelength'), ylabel('spectra')
title('Wavelength vs spectra')
subplot(2,2,2)
plot(x,x2);
axis tight
xlabel('wavelength'), ylabel('spectra')
title('Wavelength vs spectra')
subplot(2,1,2)
plot(x,x3);
axis tight
xlabel('wavelength'), ylabel('spectra')
title('Wavelength vs spectra')

 
 figure
 [t,p,r2]=pca(x1);
 plot(r2)
 figure
 [t,p,r2]=pca(x2);
 plot(r2)
 figure
 [t,p,r2]=pca(x3);
 plot(r2)
 
 %% Task 3
 % fitting sigmoid function
 clc
 load data_5_3.mat

 x=xy(:,1); y=xy(:,2);
 x_range = linspace(x(1),x(end),100);

 fun = @(tet,x) 1./(1+exp(-(tet(1)+tet(2)*x))); % Defining the exponential model

optimize_fun = @(tet)sum((y-fun(tet,x)).^2); % optimization function 
x0 = [0,0];
options = optimset('Display','iter','TolFun',1e-4,'TolX',1e-4,'MaxFunEvals',1000);
fit3 = fminsearch(optimize_fun,x0,options);
plot(x_range,fun(fit3,x_range), 'r-','LineWidth',2);
grid on;
axis tight;

 

%% Task 4
% fitting polynomial 
clc
load data_5_4.mat

x1 = x(:,1);
x2 = x(:,2);
y = yr;
x_range = linspace(x(1),x(end),9);
fun_1 = @(tet,x) (x1+x2.^2); % Defining the exponential model
optimize_fun_1 = @(tet)sum((y-fun_1(tet,x)).^2); % optimization function 
x0 = [1,1];
options = optimset('Display','iter','TolFun',1e-4,'TolX',1e-4,'MaxFunEvals',1000);
fit1 = fminsearch(optimize_fun_1,x0,options);
plot(x_range,fun(fit1,x_range), 'r-','LineWidth',2);
grid on;
axis tight;

fun_2 = @(tet,x) (x1.^2+x2); % Defining the exponential model
optimize_fun_2 = @(tet)sum((y-fun_2(tet,x)).^2); % optimization function 
x0 = [1,1];
options = optimset('Display','iter','TolFun',1e-4,'TolX',1e-4,'MaxFunEvals',1000);
fit2 = fminsearch(optimize_fun_2,x0,options);
plot(x_range,fun(fit2,x_range), 'r-','LineWidth',2);
grid on;
axis tight;



%b1 = beta(1);
%b2 = beta(2);
%beta0 = ones(3,1);
%modelfun = @(b,x) (b1*x1+b2*x1.^2);
%mdl = fitnlm(x,y,modelfun,beta0);


%% Task 5
%