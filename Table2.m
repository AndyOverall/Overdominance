% Table2.m
% Determines T1/2 for infinite and finite populations

clear all
close all
clc

% Parameters
Nset=[5e2,1e3,2e3];
u=1e-8

% Later time T12 calculations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hi=-0.01;
hf=0;

wi=1-hi;
wf=1-hf;

% Infinite population
Xi=((2*(1-hi)*u)/(hi+(2-3*hi)*u+sqrt(hi^2*(1+u)^2+4*(1-2*hi)*u)));
Xf=((2*(1-hf)*u)/(hf+(2-3*hf)*u+sqrt(hf^2*(1+u)^2+4*(1-2*hf)*u)));
A=(((1-2*hf))/((1-hf)-(1-2*hf)*Xf));
B=((1+(1-2*hf)*u+(1-2*hf)*(1-u)*Xf)/((1-u)*((1-hf)-(1-2*hf)*Xf)));
R=-(1-B)/A;
T=log((1+2*R/(Xi-Xf))/(1+R/(Xi-Xf)))/log(B)

% Finite population
for i=3:-1:1
      N=Nset(i)
      xvec=(0:N)'/2/N;
      % Calculate transition matrix Wi
      Wi=zeros(N+1,N+1);
      w=wi;
      for n=0:N
            xn=n/2/N;
            Fn=(u*w+((w-1)-u*(3*w-1))*xn-(1-u)*(2*w-1)*xn^2)/((1+u*(2*w-1))+(1-u)*(2*w-1)*xn);
            an=log(2*xn+2*Fn);
            bn=log(1-2*xn-2*Fn);
            for m=0:N
                  Wi(m+1,n+1)=exp(gammaln(N+1)-gammaln(N-m+1)-gammaln(m+1)+m*an+(N-m)*bn);
            end
      end
      
      % Determine stationary distribution
      [A,B]=eig(Wi);
      B=diag(real(B));
      [ignore,bb]=sort(B);
      b=bb(end);
      lambda=B(b);
      Phi_i=A(:,b)/sum(A(:,b));
      % Mean initial frequency in a finite population, at stationarity
      EXi=xvec'*Phi_i;
      
      % Calculate transition matrix Wf
      Wf=zeros(N+1,N+1);
      w=wf;
      for n=0:N
            xn=n/2/N;
            Fn=(u*w+((w-1)-u*(3*w-1))*xn-(1-u)*(2*w-1)*xn^2)/((1+u*(2*w-1))+(1-u)*(2*w-1)*xn);
            an=log(2*xn+2*Fn);
            bn=log(1-2*xn-2*Fn);
            for m=0:N
                  Wf(m+1,n+1)=exp(gammaln(N+1)-gammaln(N-m+1)-gammaln(m+1)+m*an+(N-m)*bn);
            end
      end
      
      % Determine stationary distribution corresponding to hf
      [A,B]=eig(Wf);
      B=diag(real(B));
      [ignore,bb]=sort(B);
      b=bb(end);
      lambda=B(b);
      Phi_f=A(:,b)/sum(A(:,b));
      % Mean final frequency in a finite population, at stationarity
      EXf=xvec'*Phi_f;
      
      % Evolve distribution
      t=0;
      EX=1;       % for later calculation
      Phi=Phi_i;
      while EX>(EXi+EXf)/2    % for later calculation
            t=t+1;
            EX=xvec'*Phi;
            Phi=Wf*Phi/sum(Wf*Phi);
      end
      T=t
end