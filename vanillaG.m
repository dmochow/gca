function G = vanillaG(x,y,L)
%UNTITLED2 granger causality from x to y with lag L
%   Detailed explanation goes here

yp = wienerAperture(y,L);
hwr = myWiener(y,yp);
mmse_r_emp = mmse(hwr,y,yp);

xp=wienerAperture(x,L);
xpyp=cat(2,xp,yp); % combine aperture of causee and causer
hwf=myWiener(y,xpyp);
mmse_f_emp = mmse(hwf,y,xpyp);

G = 1 - mmse_f_emp/mmse_r_emp;
end

