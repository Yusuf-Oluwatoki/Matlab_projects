clear all
clc
close all
addpath(genpath('./mcmcstat-master'))
%rng(0);
%%
k=6;m=2^k;n=2^k;l=1.3; nu=1.8; var=1;
s=20;
%%
[j,i] = meshgrid(1:s,1:s);
x=i(:); y=j(:); 

dmat=zeros(length(x),length(y));
for ii=1:length(x)
    for jj=1:length(y)
        deltax=(x(ii)-x(jj))/s; deltay=(y(ii)-y(jj))/s;
        dmat(ii,jj)=sqrt(deltax.^2+deltay.^2);
    end
end

%%
% d=distmat(a);
% dmat(boolean(eye(s)))=0.02;
% figure
% imagesc(dmat)
% colormap(jet)
%%
bgra=sqrt(2*nu).*(dmat./l);
Cpr= var*(2^(1-nu)/gamma(nu))*(bgra).^nu.*besselk(nu,bgra);
 %Cpr=Cpr+0.02*eye(size(Cpr));
 index = isnan(Cpr);
Cpr(index)=var;
%Cpr=Cpr+0.00001*eye(size(Cpr));
 %%
% Cpr(boolean(eye(n)))=var;
% figure
% imagesc(Cpr)
% colormap(jet)
%%
ch=chol(Cpr,'lower');
%W = randn(s*s,1);
W = randn(s*s,1);
matern_vec = ch*W;
img = reshape(matern_vec, [s, s]);

%%
figure
imagesc(img)
% colormap(jet)
colormap(jet)
%%
% n=2;
% deter= ett(Cpr);
% ett(Cpr)
%%
% params= logp(n,matern_vec,Cpr)
% logposterior=@(p) -nloglike(yfit3,p)-(p(1)-0.03)^2/(2*0.05^2)-(p(2)-0.03)^2/(2*0.05^2);
% maximum_li=fminsearch(@(p) nloglike(yfit3,p),[1,1]')
logposterior=@(p) logp(p,matern_vec,dmat)-(p(1)-0.03)^2/(2*0.05^2)-(p(2)-0.03)^2/(2*0.05^2);
maximum_li=fminsearch(@(p) logp(p,matern_vec,dmat),[1,1]')
%%
% clc
% c2=mcmc(logposterior,maximum_li,1500)
%%
ssfun= @(theta, data) logp(theta,data,dmat); %defining ssfun
params= { 
    {'b_0', maximum_li(1),0,Inf}
    {'b_1', maximum_li(2),0,Inf}
    };

model.ssfun= ssfun; data=matern_vec; model.N=1;
RR{1}= eye(2)*0.0001; RR{2}= eye(2)*0.001;
RR{3}= eye(2)*0.001; RR{4}= eye(2)*0.01;
% RR{1}= eye(2)*1.0; RR{2}= eye(2)*0.001;
% RR{3}= eye(2)*0.0001; RR{4}= eye(2)*0.0001;
options.RDR= RR; options.ntry=4;
options.nsimu=20000; options.method='ram';
[result,chain,s2chain] = mcmcrun(model,data,params,options);
%%

% clc
% c1=mcmc2(logposterior,maximum_li,500);

% %%
% figure
% subplot(2,1,1);
% plot(c2(200:end,1));hold on;
% subplot(2,1,2); 
% plot(c2(200:end,2))
% 
% figure
% plot(c2(200:end,2),c2(200:end,1))
% %%
figure
subplot(2,1,1);
plot(chain(20:end,1));hold on;
subplot(2,1,2); 
plot(chain(20:end,2))

%%
Chain=chain;

chain=Chain(5000:end,:);
% %Plot of the chains
figure
subplot(2,1,1)
plot(chain(:,1), 'Color',[0.65 0.65 0.65])
hold on
yline(mean(chain(:,1)),'k','LineWidth',3);
ylabel('$l$' , 'interpreter', 'latex')
axis tight
%grid on;
ax = gca; % current axes
ax.FontSize = 16;

subplot(2,1,2)
plot(chain(:,2), 'Color',[0.65 0.65 0.65])
hold on
yline(mean(chain(:,2)),'k','LineWidth',3);
axis tight;
%grid on;
%title('(b)')
xlabel('Length of chain' , 'interpreter', 'latex')
ylabel('$\nu$' , 'interpreter', 'latex')
x0=10;
y0=10;
width=750;
height=400;
set(gcf,'position',[x0,y0,width,height])
ax = gca; % current axes
ax.FontSize = 16;
%%

figure
subplot(2,2,1)
histogram(chain(:,1),'FaceColor',[0 0 0])
ylabel('$l$' , 'interpreter', 'latex')
hold on
xline(l,'r', 'Linewidth', 3)
%grid on;
x0=10;
y0=10;
width=350;
height=350;
set(gcf,'position',[x0,y0,width,height])
axis tight
%axis square
ax = gca; % current axes
ax.FontSize = 12;

%figure
subplot(2,2,4)
histogram(chain(:,2), 'FaceColor',[0 0 0])
ylabel('$\nu$' , 'interpreter', 'latex')
hold on
xline(nu,'r', 'Linewidth', 3)
%grid on;
axis tight
view([90 -90])

x0=10;
y0=10;
width=350;
height=350;

set(gcf,'position',[x0,y0,width,height])

%axis square
ax = gca; % current axes
ax.FontSize = 12;

%

%figure
subplot(2,2,3)
plot(chain(:,1), chain(:,2),'.', 'color', [0.65 0.65 0.65])
hold on
plot(l,nu, 'or','LineWidth',2, 'MarkerSize', 7)
hold on
plot(mean(chain(:,1)),mean(chain(:,2)),  'ok','LineWidth',2, 'MarkerSize', 7)
xlabel('$l$' , 'interpreter', 'latex')
ylabel('$\nu$' , 'interpreter', 'latex')
legend('Chain pairs','True parameter','Estimated parameter' , 'interpreter', 'latex')
%grid on;
axis tight
%axis square
x0=10;
y0=10;
width=350;
height=350;
set(gcf,'position',[x0,y0,width,height])

ax = gca; % current axes
ax.FontSize = 12;

%%


figure
subplot(2,1,1)
[w,x,y]=autocorr(chain(:,1),'NumLags', 399,'NumStd',2);
plot(x,w,'k')
hold on
yline(y(1),'r')
hold on 
yline(y(2),'r')
grid on
ylabel('ACF $l$' , 'interpreter', 'latex')
ax = gca; % current axes
ax.FontSize = 16;
subplot(2,1,2)
[w,x,y]=autocorr(chain(:,2),'NumLags', 399,'NumStd',2);
plot(x,w,'k')
hold on
yline(y(1),'r')
hold on 
yline(y(2),'r')
grid on
xlabel('Lag' , 'interpreter', 'latex')
ylabel('ACF $\nu$' , 'interpreter', 'latex')
x0=10; y0=10; width=750; height=400;
set(gcf,'position',[x0,y0,width,height]);
ax = gca; ax.FontSize = 16;

%%
chainstats(chain,result)
%%
% clc
% x0=[2,3]; N=100; vari=0.3;
% Normal MH
 function chain=mcmc(logposterior,x0,N)
compo=numel(x0);
chain=zeros(N,compo);
chain(1,:)=x0;
x=x0;
acc=0;
val=logposterior(x0);
R=eye(compo)*0.0001;
R=chol(R);
for i=2:N
    proposal=x+R*randn(compo,1);
    newval=logposterior(proposal);
    if (log(rand()) < val-newval)
        x=proposal;
        val=newval;
        acc=acc+1;
    end
    chain(i,:)=x;
end
ratio=acc/(N-1)
 end
 
 % Adaptive MH
 function chain=mcmc2(logposterior,x0,N)
compo=numel(x0);
chain=zeros(N,compo);
chain(1,:)=x0; n0=1;
x=x0; t0=20; s0=8;
acc=0;
val=logposterior(x0);
R=eye(compo)*0.0001;
R=chol(R);
for i=2:N
    proposal=x+R*randn(compo,1);
    newval=logposterior(proposal);
    if (log(rand()) < val-newval)
        x=proposal;
        val=newval;
        acc=acc+1;
    end
    chain(i,:)=x;
    %is2=gamrnd((n0+10)/2, (n0*s0^2+sumsqr(x))/2);
    %s2(i)=1/is2;
    if (i>t0)
        C=cov(chain(1:i,:))+eye(compo)*0.000001;
        %C=2.4^2/compo*(C+eye(compo)*0.000001*(abs(max(C(:))) < 0.00001));
        R=chol(C,'lower');
    end     
end
ratio=acc/(N-1)
 end
%%
function params= logp(p,img,dmat)
var=1;
bgra=sqrt(2*p(2)).*(dmat./p(1));
Cpr= var*(2^(1-p(2))/gamma(p(2)))*(bgra).^p(2).*besselk(p(2),bgra);
index = isnan(Cpr);
Cpr(index)=var+1e-07;
spar=Cpr\speye(size(Cpr));
params=-0.5.*(img'*spar*img) +(1/2*log(2*pi) -0.5.*ett(Cpr));
params=-params;

end

function deter= ett(Cpr)
        chhe=chol(Cpr, 'lower');
        deter= 2*sum(log(diag(chhe)));       
end
function likelihood=nloglike(x,p)
likelihood=0;

for i=1:numel(x)-1
    likelihood=likelihood + 1/2*log(2*pi*p(1)/(2*p(2))...
        *(1-exp(-2*p(2))))+ p(2)/p(1)*(1-exp(-2*p(2)))*(x(i+1)-exp(-p(2))*x(i))^2;
end
end