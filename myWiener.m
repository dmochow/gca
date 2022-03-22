function [hw,r,R] = myWiener(xd,Xp)
%[HW,r,R]=MYWIENER(xd,Xp) PREDICT xD FROM APERTURE Xp


% xd = desired signal (time series with dim 1)
% Xp = prediction aperture (multivariate time series)
xd=xd(:);
N=numel(xd);
L=size(Xp,2); % this is not necessarily the max lag L but we still call this variable L here
r = mean ( repmat(xd,1,L) .* Xp );  % prediction aperture
r=r(:); % force column vector
R = 1/N*(Xp.'*Xp); % temporal covariance matrix of prediction aperture
hw=R\r; % Wiener filter
end

