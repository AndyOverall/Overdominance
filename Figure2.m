clc
close all

load Data1

for mu=1:2
      u=uset(mu);
      Xhat=sqrt(u)/(1+sqrt(u));
      Xt(T1,mu)=Xhat;
      x=Xhat;
      
      for t=T2
            h=hstar;
            F=(u*(1-h)-(h-u*(2-3*h))*x-(1-u)*(1-2*h)*x^2)/((1+u*(1-2*h))+(1-u)*(1-2*h)*x);
            Xt(t,mu)=x+F;
            x=Xt(t,mu);
      end
      
      for t=T3
            h=0;
            F=(u*(1-h)-(h-u*(2-3*h))*x-(1-u)*(1-2*h)*x^2)/((1+u*(1-2*h))+(1-u)*(1-2*h)*x);
            Xt(t,mu)=x+F;
            x=Xt(t,mu);
      end
end

hold on

lw=3;
mu=2;
plot(T,log10(Xt(:,mu)),  'k','linewidth',lw)
plot(T,log10(EX(:,mu,3)),'r','linewidth',lw)
plot(T,log10(EX(:,mu,2)),'b','linewidth',lw)
plot(T,log10(EX(:,mu,1)),'m','linewidth',lw)

mu=1;
thin=50:170:Lt;
lw=3.01;

plot(T(thin),log10(Xt(thin,mu))  ,'k-o','linewidth',lw)
plot(T(thin),log10(EX(thin,mu,3)),'r-o','linewidth',lw)
plot(T(thin),log10(EX(thin,mu,2)),'b-o','linewidth',lw)
plot(T(thin),log10(EX(thin,mu,1)),'m-o','linewidth',lw)

v=(2:-0.01:-7);
vv=ones(size(v));
plot(0*vv,v,'k','linewidth',1)
plot(2000*vv,v,'k','linewidth',1)

hl=legend(' {\itu} = 10^{-5}, {\itN} = \infty{ }',...
          ' {\itu} = 10^{-5}, {\itN} = 2000{ }',...
          ' {\itu} = 10^{-5}, {\itN} = 1000',...
          ' {\itu} = 10^{-5}, {\itN} = 500',...
          ' {\itu} = 10^{-8}, {\itN} = \infty',...
          ' {\itu} = 10^{-8}, {\itN} = 2000',...
          ' {\itu} = 10^{-8}, {\itN} = 1000',...
          ' {\itu} = 10^{-8}, {\itN} = 500','location','eastoutside');

f=25;
h=-4.55;
text(-700-250,h,'{\ith} = 0','fontsize',f)
text(1000-830,h,'{\ith} = - 0.01','fontsize',f)
text(3000-460,h,'{\ith} = 0','fontsize',f)

axis([-1000,4000,-6.5,-1.7])
set(gca,'linewidth',3,'fontsize',20)
xlabel('time, \itt','fontsize',25)
ylabel('log_{10}({\itE }[{\itX_t} ])','fontsize',25)
set(gca,'linewidth',3,'fontsize',25)
set(gca,'xtick',[0,2000,4000],'ytick',[-6:-2])
orient('landscape')
box on