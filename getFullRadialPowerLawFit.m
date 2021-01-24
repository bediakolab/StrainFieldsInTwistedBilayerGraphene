function [rms_deviation] = getFullRadialPowerLawFit(c,fitdata)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
matsize = size(fitdata);
predmat = getFullRadialPowerLawPred(c,matsize);
diff = predmat - fitdata;
diff(fitdata == -1) = 0;  % because -1 was the signal that this point should be excluded.
% Compute rms manually, correcting for the number of excluded points.
num_excluded = nnz(fitdata == -1);
num_total = numel(fitdata);
num_to_use = num_total - num_excluded;
rms_deviation = sqrt(sum(sum(diff.^2))/num_to_use);

end

