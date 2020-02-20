% Figure1

clear all
close all
clc

% Big plot
h=-1:0.001:1;
w=1-h;
u=1e-4;
Xhat=2*u*w./((1-w)-u*(1-3*w)+sqrt((1+u^2)*(1-w).^2+2*u*(w.^2+2*w-1)));
plot(h,Xhat,'b','linewidth',3)
hold on
u=1e-8;
Xhat=2*u*w./((1-w)-u*(1-3*w)+sqrt((1+u^2)*(1-w).^2+2*u*(w.^2+2*w-1)));
plot(h,Xhat,'r','linewidth',3)
hold on
axis([-1,1,-0.01,0.4])
hl=legend(' {\itu} = 10^{-4},{ }',' {\itu} = 10^{-8}{ }',...
          'location','north','Orientation','horizontal');
set(hl,'fontsize',30)

u=1e-4;
Xhat=2*u*w./((1-w)-u*(1-3*w)+sqrt((1+u^2)*(1-w).^2+2*u*(w.^2+2*w-1)));
plot(h,Xhat,'b','linewidth',3)

xlabel('\ith','fontsize',30)
ylabel('$$\hat{X}$$','Interpreter','Latex','fontsize',30)
set(gca,'linewidth',3,'fontsize',20)
set(gca,'ytick',[0,0.1,0.2,0.3],'yticklabel',{'0','0.1','0.2','0.3'})
set(gca,'xtick',[-1,-0.5,0,0.5, 1],'xticklabel',{'-1','-0.5','0','0.5','1'})

orient('landscape')

% Small plot
h=-0.02:0.0001:0.02;
w=1-h;
axes('Position',[0.575,0.26,0.27,0.5])

box on

u=1e-4;
Xhat=2*u*w./((1-w)-u*(1-3*w)+sqrt((1+u^2)*(1-w).^2+2*u*(w.^2+2*w-1)));
plot(h,Xhat,'b','linewidth',3)
hold on
u=1e-8;
Xhat=2*u*w./((1-w)-u*(1-3*w)+sqrt((1+u^2)*(1-w).^2+2*u*(w.^2+2*w-1)));
plot(h,Xhat,'r','linewidth',3)
hold on

axis([-0.02,0.02,-0.0005,0.022])
set(gca,'ytick',[0,0.01,0.02],'yticklabel',{'0','0.01','0.02'})
set(gca,'xtick',[-0.02,0,0.02],'xticklabel',{'-0.02','0','0.02'})


set(gca,'linewidth',3,'fontsize',18)

orient('landscape')
ounits=get(gcf,'Units');
set(gcf,'Units','normalized','OuterPosition',[0 0 1 1],'Units',ounits)