% Table1.m
% Determines equilbrium and expected properties of the mean allele frequency

clear all
close all
clc

% Parameters
Nset=[5e2,1e3,2e3];
N=Nset(1);
u=1e-5;
hset=[-0.1,-0.01,0,0.01,0.1]';

%Useful
EX=zeros(5,1);
Xhat=zeros(5,1);

for i=1:5
      i
      h=hset(i);
      w=1-h;
      % Calculate transition matrix
      W=zeros(N+1,N+1);
      for n=0:N
            xn=n/2/N;
            Fn=(u*w+((w-1)-u*(3*w-1))*xn-(1-u)*(2*w-1)*xn^2)/((1+u*(2*w-1))+(1-u)*(2*w-1)*xn);
            an=log(2*xn+2*Fn);
            bn=log(1-2*xn-2*Fn);
            for m=0:N
                  W(m+1,n+1)=exp(gammaln(N+1)-gammaln(N-m+1)-gammaln(m+1)+m*an+(N-m)*bn);
            end
      end
      
      % Determine stationary distribution
      [A,B]=eig(W);
      B=diag(real(B));
      [a,bb]=sort(B);
      b=bb(end);
      lambda=B(b)
      phi=A(:,b)/sum(A(:,b));
      xvec=(0:N)'/2/N;
      
      % Calculate infinite population equilibrium frequency
      Xhat(i)=2*u*w/((1-w)-u*(1-3*w)+sqrt((1+u^2)*(1-w)^2+2*u*(w^2+2*w-1)));
      
      % Calculate mean frequency in a finite population, at stationarity
      EX(i)=xvec'*phi;
end
u
N
log10(EX)
log10(Xhat)


