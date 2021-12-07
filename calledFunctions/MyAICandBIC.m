function [AIC,BIC] = MyAICandBIC(MSE,n,p)

%n=ndata
%MSE=mean squared error
%p=number of paramters

MSE=MSE(:);
n=n(:);
p=p(:);

AIC=log(MSE).*n + p.*2 + (p.*(p+1).*2)./(n-p-1);
BIC=log(MSE).*n + log(n).*p;
