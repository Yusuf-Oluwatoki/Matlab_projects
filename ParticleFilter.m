%% Particle filter
clear all
close all
clc
rng('default')
%% add toolbox path and load mcmc chain for comparison
load('eam.mat') % AM with ode45 (can be changed to eam, edram, dram, mh, emh
               % eram) e means chains from Euler method.
addpath(genpath('./mcmcstat-master'))
%% Parameters for the forced spring-mass system

p=[4 1]; %v is 2 and gamma is 1
s=4;    %  ode solver to be used (1=ode23,2=ode45,3=ode113,4=Euler)
h=0.05; T=15; t=0:h:T; %time steps

x=[1;0]; y0 = [1;0.001];  %initial values for noisy system and true system
N=length(t); %length of the system

[A,state]=forced_spring_ode(t,x,p,s); %solving the system with ode solvers

n=2; nobs=1; %numbers of states and number of observation per time

M=[0 1; -p(1) -p(2)]; %state matrix coefficient
A= eye(n)+ M*h;       %state matrix
H=[1 0]; L=h*[0 0; 0 1]; %measurement matrix
Q=0.03^2*[1 0; 0 1]; % covariance of state
q=chol(Q, 'lower'); %std of state
R= 0.05^2*eye(nobs);  % covariance of measurement
r=chol(R,'lower'); %std of measurement
state=state'; %transpose state to row n by N matrix

IR=inv(R);

f=@(x) [x(1)+ h*x(2); -p(1)*x(1)*h-(1+p(2))*x(2)*h];
h=@(x) x(1);

y=zeros(nobs,N); %initialize measurement
y(:,1)= H*state(:,1)+r*randn(nobs,1); %first value of noisy measurement


for i=2:N
     y(:,i)=h(state(:,i-1))+r*randn(nobs,1); %noisy measurements
end 

%% Particle filter steps
m=zeros(n,N); 
m(:,1)=[1,0];
Npart=700;
particle=5*randn(n,Npart);


for i=2:N
   w= zeros(Npart,1);
   for j=1:Npart
       particle(:,j)=f(particle(:,j))+q*randn(n,1);
       w(j)=exp(-0.5*(h(particle(:,j))-y(:,i))'*IR*(h(particle(:,j))-y(:,i)));
   end
   wtot=sum(w);
   w=w./wtot;
   
   ind=randsample(1:Npart,Npart,true,w);
   particle=particle(:,ind);
   m(:,i)= mean(particle,2);
end
%% Root mean squared error of kalman filter estimates and state values
rmse=sqrt(mean((m(1,N:end)-y(1,N:end)).^2)) %paticle filter and noisy system

%MCMC mean used in the ode solver to compare with kalman filter estimates
[~,yfit2]= forced_spring_ode(t,y0,[mean(chain(:,1)) mean(chain(:,2))],s);
yfit2=yfit2';
rmse2=sqrt(mean((m(1,N:end)-yfit2(1,N:end)).^2)) %kalman filter and MCMC mean

rmse3=sqrt(mean((y(1,N:end)-yfit2(1,N:end)).^2)) %noisy system and MCMC mean


%% Plots of estimates and system fit
figure
hold on   
%plot(t,state(1,:),'k', 'linewidth',1)
hold on
plot(t,y(1,:),'k--', 'linewidth',1)
hold on
plot(t,m(1,:),'r--', 'linewidth',1)
hold on 
plot(t,yfit2(1,:),'k', 'linewidth',1)
legend('Noisy measurement','Particle filtering estimate','MCMC mean','interpreter','latex')

xlabel('Time', 'interpreter', 'latex', 'fontsize', 20);
ylabel('State', 'interpreter', 'latex', 'fontsize', 20);
title('Forced spring-mass system', 'interpreter', 'latex', 'fontsize', 20)
box on
axis tight
grid on
x0=10;
y0=10;
wihh=600;
height=400;
set(gcf,'position',[x0,y0,wihh,height]);
ax=gca;
ax.FontSize = 20;


% 
% rmse=sqrt(mean(state(1,floor(N/2):end)-m(1,floor(N/2):end)).^2);
% 
% [~,yfit2]= forced_spring_ode(t,y0,[mean(chain(:,1)) mean(chain(:,2))],s);
% 
% %%
% figure
% hold on   
% %plot(t,state(1,:),'k', 'linewidth',1)
% hold on
% plot(t,y(1,:),'k--', 'linewidth',1)
% hold on
% plot(t,m(1,:),'r--', 'linewidth',1)
% hold on 
% plot(t,yfit2(:,1),'k', 'linewidth',1)
% legend('Noisy measurement','Particle filtering estimate','MCMC mean','interpreter','latex')
% 
% xlabel('Time', 'interpreter', 'latex', 'fontsize', 20);
% ylabel('State', 'interpreter', 'latex', 'fontsize', 20);
% title('Noisy forced spring-mass system', 'interpreter', 'latex', 'fontsize', 20)
% box on
% axis tight
% grid on
% x0=10;
% y0=10;
% width=600;
% height=400;
% set(gcf,'position',[x0,y0,width,height])
% ax=gca;
% ax.FontSize = 20;
% 
% 

