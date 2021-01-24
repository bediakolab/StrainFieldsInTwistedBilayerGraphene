function [fit1,optimalparams1] = getPowerLawFit(ydata,xdata,coeff_guesses)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

disp('Select the left x of the background region 1.');
[x1,y1] = ginput(1);
disp('Select the right x of the background region 1.');
[x2,y2] = ginput(1);

fitdata1 = [xdata',ydata'];
fitdata1(fitdata1(:,1) < x1,:) = [];
fitdata1(fitdata1(:,1) > x2,:) = [];


powerlaw_pred = @(c,x) c(1)*(x).^(c(2));
powerlaw_fit1 = @(c) rms(fitdata1(:,2) - powerlaw_pred(c,fitdata1(:,1)));
options = optimset;
options.MaxFunEvals = 10000;
optimalparams1 = fminsearch(powerlaw_fit1,coeff_guesses);
fit1 = powerlaw_pred([optimalparams1(1),optimalparams1(2)],xdata);

end

