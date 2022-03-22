function [mmse,xhat] = mmse(hw,x,xp)
%MMSE Wiener filter minimum mean squared error when predicting time series 
% x from prediction vector xp using temporal filter hw
% 
x=x(:);
hw=hw(:);
if numel(hw)~=size(xp,2), error('number of elements in hw must match number of columns in xp'); end
if size(x,1)~=size(xp,1), error('x and xp must have same # of rows'); end
xhat=xp*hw; 
mmse=mean((x-xhat).^2);
end

