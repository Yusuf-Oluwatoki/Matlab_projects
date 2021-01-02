clear all
close all
clc

%% Add toolbox path
addpath(genpath('./mcmcstat-master'))

%% Parameters for the forced spring-mass system

p=[4 1]; %v is 2 and gamma is 1

x=[1;0]; % initial point for solution
h=0.1; T=15; t=0:h:T; %time steps

[A,state]=forced_spring_ode(t,x,p,1); % system solution with ode23 s=1
[A2,state2]=forced_spring_ode(t,x,p,2);% system solution with ode45 s=2
[A3,state3]=forced_spring_ode(t,x,p,3);% system solution with ode113 s=3
[tt,state4]=eul(t,x,p);                 % system solution with Euler s=4 or nothing

obnoise=0.05;  %observation noise

obs=state(:,1)+obnoise*randn(size(state(:,1))); % noisy system with ode23
obs2=state2(:,1)+obnoise*randn(size(state2(:,1))); % noisy system with ode45
obs3=state3(:,1)+obnoise*randn(size(state3(:,1))); % noisy system with ode113
obs4=state4(:,1)+obnoise*randn(size(state4(:,1)));% noisy system with Euler



%% Plot of all state solutions
figure; hold on;
plot(t,state(:,1),'k*','linewidth',1);
plot(t,state2(:,1),'k','linewidth',1);
plot(t,state3(:,1),'k--','linewidth',1);
plot(t,state4(:,1),'ko','linewidth',1);
legend('Ode23','Ode45','Ode113','Euler', 'interpreter', 'latex', 'fontsize', 20)
xlabel('Time', 'interpreter', 'latex', 'fontsize', 20);
ylabel('State', 'interpreter', 'latex', 'fontsize', 20);
title('Forced spring-mass system', 'interpreter', 'latex', 'fontsize', 20)
box on
axis tight
grid on
x0=10;
y0=10;
width=600;
height=400;
set(gcf,'position',[x0,y0,width,height])
ax=gca;
ax.FontSize = 20;

%% Plots of noisy systems
figure
sgtitle('Noisy forced spring-mass system', 'interpreter', 'latex', 'fontsize', 20)
ax=gca;
ax.FontSize = 20;
subplot(2,2,1); plot(t,obs,'k','linewidth',1);grid on;axis tight; legend('Ode23','interpreter', 'latex', 'fontsize', 20)
subplot(2,2,2); plot(t,obs2,'r','linewidth',1);grid on;legend('Ode45','interpreter', 'latex', 'fontsize', 20)
subplot(2,2,3); plot(t,obs3,'k','linewidth',1);grid on; legend('Ode113','interpreter', 'latex', 'fontsize', 20)
subplot(2,2,4); plot(t,obs4,'r','linewidth',1);grid on; legend('Euler','interpreter', 'latex', 'fontsize', 20)
axis tight
x0=10;
y0=10;
width=600;
height=400;
set(gcf,'position',[x0,y0,width,height])