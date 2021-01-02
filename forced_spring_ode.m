function [A,state]=forced_spring_ode(time,xinit,p,s)
% ode solver for the forced spring-mass system.
% Function takes in the spring system time, initial values, parameters v and gamma
% as a vector p and s as the soltion method.
%
% 

k=p(1); b=p(2); %assign parameter to the variables in the function

A=[0 1;-k -b];  %system matrix coefficients
rhs=@(t,x) A*x; %system matrix form


%% Make choice for ode solvers
if s==1  %ode23 for s==1
    [~,state]=ode23(rhs,time,xinit);
elseif (s==2) %ode45 for s==2
    [~,state]=ode45(rhs,time,xinit);
elseif (s==3) %ode113 for s==3
    [~,state]=ode113(rhs,time,xinit);
else %Euler method for s any number or nothing
    [~,state]=eul(time,xinit,p);
end
end