%% MCMC estimation
clear all
close all
clc
%% Add toolbox path
addpath(genpath('./mcmcstat-master'))
%% Parameters for the forced spring-mass system

p=[4 1]; %v is 2 and gamma is 1
s=4;    %  ode solver to be used (1=ode23,2=ode45,3=ode113,4=Euler)
h=0.05; T=15; t=0:h:T; %time steps

x=[1;0]; y0 = [1;0.001];  %initial values for noisy system and true system
N=length(t); %length of the system

[A,state]=forced_spring_ode(t,x,p,s); %solving the system with ode solvers

obnoise=0.05;  %observation noise

obs=state(:,1)+obnoise*randn(size(state(:,1))); % noisy system with any ode solver

%% Plot of system and noise
figure
plot(t,obs,'k--','markerfacecolor','b','markersize',1)
hold on
plot(t,state(:,1),'*k','linewidth',1);
xlabel('Time','fontsize',20, 'interpreter','latex');
ylabel('State','fontsize',20, 'interpreter','latex');
legend('Noisy system', 'True system', 'interpreter','latex')
title('Forced spring-mass system', 'interpreter','latex')
box on
axis tight
grid on;
x0=10;
y0=10;
width=600;
height=400;
set(gcf,'position',[x0,y0,width,height])
ax=gca;
ax.FontSize = 20;
%% State optimization from noisy data

%put noisy data in structure format
state_noise.ydata(:,1)=t; state_noise.ydata(:,2)=obs;
state_noise.y0 = [1;0.001]; %initial value for noise data

k00 = [10 3]'; %initial guess of parameter
[fit,rss] = fminsearch(@forced_spring_ss,k00,[],state_noise,s) %optimal parameter search
mse = rss/(length(state_noise.ydata)-2); %mean squared error of previous system and new system

[~,yfit]= forced_spring_ode(t,state_noise.y0,fit,s); %new system values

%% MCMC runs using toolbox
%initialising runs with optimal parameters searched above 
params = {
    {'k1', fit(1), 0}
    {'k2', fit(2), 0}
    };
data=state_noise; %noisy system to be estimated

ssfun= @(theta,data,s)forced_spring_ss(theta,data,s); %defining ssfun
model.ssfun = ssfun;
model.sigma2 = mse;
options.updatesigma = 1;
nsimu=10000;  %chain length
options.nsimu=nsimu; options.method='am'; %adaptive method
[result,Chain,s2chain] = mcmcrun(model,data,params,options); %chain run

chain=Chain(2000:end,:); %remove burn-in from chain
chainstats(chain,result) %chain statistics

%% Predictive envelopes from chain

yfit3=zeros(nsimu,length(t));
for i=1:length(chain)
    [~,state_all]= forced_spring_ode(t,state_noise.y0,chain(i,:),s);
    yfit3(i,:)=state_all(:,1);
end
%chain mean estimate
[~,yfit2]= forced_spring_ode(t,state_noise.y0,[mean(chain(:,1)) mean(chain(:,2))],s);

%% Plot of predictive envelope and chain mean estimates
figure
plot(t,state_noise.ydata(:,2),'k--')
hold on
plot(t,yfit2(:,1),'r-');
hold on
plot(t,yfit(:,1),'k-');
hold on
plot(t,yfit3, 'Color',[0.9 0.9 0.9])
hold on
plot(t,yfit2(:,1),'r-','linewidth',1);
hold on
plot(t,yfit(:,1),'k-','linewidth',1);
hold on
plot(t,state_noise.ydata(:,2),'k--','linewidth',1)
xlabel('Time', 'interpreter', 'latex', 'fontsize', 20)
ylabel('State', 'interpreter', 'latex', 'fontsize', 20)
legend('Noisy data', 'MCMC mean','Optimization fit','MCMC runs', 'interpreter', 'latex', 'fontsize', 20)
title('Predictive envelopes', 'interpreter', 'latex', 'fontsize', 20);
axis tight
grid on;
x0=10;
y0=10;
width=600;
height=400;
set(gcf,'position',[x0,y0,width,height])
ax=gca;
ax.FontSize = 20;


%% Plot of the chains
figure
subplot(2,1,1); plot(chain(:,1), 'Color',[0.65 0.65 0.65]);hold on
yline(mean(chain(:,1)),'k','LineWidth',3);
ylabel('$v^2$' , 'interpreter', 'latex', 'fontsize', 20);axis tight;grid on
ax=gca;
ax.FontSize = 20;
subplot(2,1,2)
plot(chain(:,2), 'Color',[0.65 0.65 0.65])
hold on
yline(mean(chain(:,2)),'k','LineWidth',3);
axis tight;
grid on;
xlabel('Length of chain' , 'interpreter', 'latex', 'fontsize', 20)
ylabel('$g$' , 'interpreter', 'latex', 'fontsize', 20)
x0=10;
y0=10;
width=600;
height=350;
set(gcf,'position',[x0,y0,width,height])
ax=gca;
ax.FontSize = 20;
%% plots of chain distribution
figure

subplot(2,2,1);
histogram(chain(:,1),'FaceColor',[0 0 0]); xline(p(1),'r', 'Linewidth', 3);
grid on;ax=gca;ax.FontSize = 20;
subplot(2,2,4);
histogram(chain(:,2), 'FaceColor',[0 0 0]);xline(p(2),'r', 'Linewidth', 3);
grid on; view([90 -90]);ax=gca; ax.FontSize = 20;
subplot(2,2,3);
plot(chain(:,1), chain(:,2),'.', 'color', [0.65 0.65 0.65]);hold on;
plot(p(1),p(2), 'or','LineWidth',2, 'MarkerSize', 7);
plot(fit(1),fit(2),  'ok','LineWidth',2, 'MarkerSize', 7);grid on
xlabel('$v^2$' , 'interpreter', 'latex')
ylabel('$g$' , 'interpreter', 'latex')
legend('Chain pairs','True parameter','Estimated parameter' , 'interpreter', 'latex','Location','best')

x0=10;
y0=10;
width=600;
height=600;
set(gcf,'position',[x0,y0,width,height])
ax=gca;
ax.FontSize = 20;
