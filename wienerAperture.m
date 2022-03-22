function xp = wienerAperture(x,L)
%R=WIENER(X,L) RETURN PREDICTION APERTURE OF TEMPORAL WIENER FILTER WITH
%MAX LAG L AND TIME SERIES X
%   

% old code, pre 2022
%xp=tplitz(x,L+1); 
%xp=xp(:,2:end-1);

%xp=tplitz1(x,L);  % uses array processing
xp=tplitz2(x,L); % uses for loop

end

