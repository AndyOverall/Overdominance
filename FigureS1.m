% FigureS1.m
% Compares two different simulations
% 1) Conditioned model, where population size is held fixed:X'=Bin(N,2X+2F)/2N
% 2) Naive extension of the weak selection model X'=Bin(2N,X+F)/2N

clear all
close all
clc

% Parameters
N=1e2;              % Census population size
u=1e-5;             % Mutation rate
h=-0.2;              % Dominance coefficient
T=2e3;              % Duration of simulation
R=1e5;              % Number of replicates
n=1;                % Initial number of disease alleles

% Useful
X1=zeros(T+1,R);
X2=zeros(T+1,R);

% Initialisation
X1(1,:)=n/2/N;
X2(1,:)=n/2/N;

% Iterate all R replicates at each time for Case 1 of strong selection
for t=1:T
    x=X1(t,:);
    F=(((1-h)*u-(h+(2-3*h)*u)*x-(1-2*h)*(1-u)*x.^2)./((1+(1-2*h)*u)+(1-2*h)*(1-u)*x));
    X1(t+1,:)=binornd(N,2*x+2*F,[1,R])/2/N;
end

% Iterate all R replicates at each time for Case 2 of weak selection
for t=1:T
    x=X2(t,:);
    F=(((1-h)*u-(h+(2-3*h)*u)*x-(1-2*h)*(1-u)*x.^2)./((1+(1-2*h)*u)+(1-2*h)*(1-u)*x));
    X2(t+1,:)=binornd(2*N,x+F,[1,R])/2/N;
end

plot((0:T),mean(X1,2),'k','linewidth',8)
hold on
plot((0:T),mean(X2,2),'k','linewidth',3)

xlabel('time, \itt','fontsize',30)
ylabel('{\itE} [{\itX_t} ]','fontsize',30)
set(gca,'linewidth',3,'fontsize',25)
set(gca,'xtick',[0,1e3,2e3],'ytick',[0,0.02,0.04,0.06])
axis([0,2e3,0,0.07])
hl=legend('  strong lethal selection{  }',...
          '  weak selection','location','south');
set(hl,'fontsize',30,'linewidth',3)
orient('landscape')
ounits=get(gcf,'Units');
set(gcf,'Units','normalized','OuterPosition',[0 0 1 1],'Units',ounits)