%Figure3Generator.m
clear all
close all
clc

%Parameters
Nset=[500,1000,2000];
uset=[1e-8,1e-5];
hstar=-0.01;

T1=[1:1000];
T2=[1001:3000];
T3=[3001:5001];
Lt=length(T1)+length(T2)+length(T3);

%Useful
EX=zeros(Lt,2,3);
Xhat=zeros(2,1);

for mu=1:2
      u=uset(mu)
      Xhat(mu)=2*(1-hstar)*u/(hstar+(2-3*hstar)*u+sqrt(hstar^2*(1+u)^2+4*(1-2*hstar)*u));
      for n=1:3
            N=Nset(n)
            xn=(0:N)'/2/N;
            
            %Transition matrix for h=0
            W0=zeros(N+1,N+1);
            Fn=(u+(-u*2)*xn-(1-u)*xn.^2)./((1+u)+(1-u)*xn);
            an=log(2*xn+2*Fn);
            bn=log(1-2*xn-2*Fn);
            for m=0:N
                  W0(m+1,:)=exp(gammaln(N+1)-gammaln(N-m+1)-gammaln(m+1)+m*an+(N-m)*bn);
            end
            
            %Transition matrix for w=1-hstar
            Wstar=zeros(N+1,N+1);
            w=1-hstar;
            Fn=(u*w+((w-1)-u*(3*w-1))*xn-(1-u)*(2*w-1)*xn.^2)./((1+u*(2*w-1))+(1-u)*(2*w-1)*xn);
            an=log(2*xn+2*Fn);
            bn=log(1-2*xn-2*Fn);
            for m=0:N
                  Wstar(m+1,:)=exp(gammaln(N+1)-gammaln(N-m+1)-gammaln(m+1)+m*an+(N-m)*bn);
            end
            
            %Stationary distribution for h=0
            [A,B]=eig(W0);
            B=diag(real(B));
            [ignore,b]=sort(B);
            b=b(end);
            lambda=B(b);
            Phi=A(:,b)/sum(A(:,b));
            EX0=xn'*Phi;
            EX(T1,mu,n)=EX0;
            
            %Evolve distribution with Wstar and determine mean
            for t=T2
                  EX(t,mu,n)=xn'*Phi;
                  Phi=Wstar*Phi/sum(Wstar*Phi);
            end
            
            %Evolve distribution with W0 and determine mean
            for t=T3
                  EX(t,mu,n)=xn'*Phi;
                  Phi=W0*Phi/sum(W0*Phi);
            end
      end
end
T=[T1,T2,T3]-1001;

save Data1