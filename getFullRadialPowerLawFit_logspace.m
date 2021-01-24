function [returnval_outer] = getFullRadialPowerLawFit_logspace(c,fitdata,optimize)
% return

% need to pre_filter the fitdata before log transform to avoid nans:
fitdata(fitdata <= 0) = 0.1;  % it's okay because they will be disregarded in the final RMSR.
% However we must remove them in the linear fit.
fitdata_log10 = log10(fitdata);

matsize = size(fitdata);
returnval = getFullRadialPowerLawPred_logspace(c,matsize,fitdata_log10,optimize);
if optimize
    predmat_log10 = returnval;
else
    predmat_log10 = returnval.predmat_log10;
end
diff = predmat_log10 - fitdata_log10;
diff(fitdata_log10 == -1) = 0;  % because -1 was the signal that this point should be excluded.
% Compute rms manually, correcting for the number of excluded points.
num_excluded = nnz(fitdata <= 0);  % 0 wrecks havoc with a log transform; eliminate it if it exists.
num_total = numel(fitdata);
num_to_use = num_total - num_excluded;
rms_deviation = sqrt(sum(sum(diff.^2))/num_to_use);

if optimize
    returnval_outer = rms_deviation;
else
    returnval_outer.RMSR = rms_deviation;
    returnval_outer.log10A = returnval.log10A;
    returnval_outer.B = returnval.B;
end

end

