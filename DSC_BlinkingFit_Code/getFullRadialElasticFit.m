function [rms_deviation] = getFullRadialElasticFit(c,fitdata,function_string,const_vals)
% Helper function for the FourDSTEM_Analysis_Engine.m class.
% Allows fitting of multiple functional forms through the string argument
% to getFullRadialElasticPred.m
% This function does not handle exceptions.
%
% Nathanael Kazmierczak, 03/12/2020
if nargin < 4
    const_vals = [];
end

matsize = size(fitdata);
predmat = getFullRadialElasticPred(c,matsize,function_string,const_vals);
diff = predmat - fitdata;
diff(fitdata == -1) = 0;  % because -1 was the signal that this point should be excluded.
% Compute rms manually, correcting for the number of excluded points.
num_excluded = nnz(fitdata == -1);
num_total = numel(fitdata);
num_to_use = num_total - num_excluded;
rms_deviation = sqrt(sum(sum(diff.^2))/num_to_use);

end

