function [rho,pval] = corrcoef12(x,y)
%[RHO,PVAL]=CORRCOEF12(X,Y)
% wrapper to grab the off-diagonal elements from corrcoef of two signals
x=x(:);
y=y(:);
[r,p]=corrcoef(x,y);
rho=r(1,2);
pval=p(1,2);
end

