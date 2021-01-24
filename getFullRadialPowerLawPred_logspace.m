function [returnval] = getFullRadialPowerLawPred_logspace(c,matsize,fitdata_log10,optimize)
% Assumes a spacing of one in each pixel direction.
% c(1) is x0, c(2) is y0.
% exponent and log10(prefactor) are calculated linearly.
% We now need fitdata in this function to do the linear bit.



x0 = c(1);
y0 = c(2);
xbase = 1:matsize(2);
ybase = 1:matsize(1);
[xspace,yspace] = meshgrid(xbase,ybase);
xspace_centered = xspace - x0;
yspace_centered = yspace - y0;
radii = sqrt(xspace_centered.^2 + yspace_centered.^2);
% Need to unwrap this to do it properly by matrices.
radii_unwrapped = radii(:);
designones = ones(size(radii));
designones_unwrapped = designones(:);
org_design = horzcat(designones_unwrapped,log10(radii_unwrapped));

fitdata_log10_unwrapped = fitdata_log10(:);
% find indices of the ones that can't be fit logarithmically.
indices = find(fitdata_log10_unwrapped == -1);
fitdata_log10_unwrapped(indices) = [];
% beause above is 1D vector, these work as subscripted indices too.
design = org_design;
design(indices,:) = [];
linear_coeffs = design\fitdata_log10_unwrapped;
predmat_unwrapped = org_design*linear_coeffs;
predmat_log10 = reshape(predmat_unwrapped,size(fitdata_log10));

if optimize
    returnval = predmat_log10;
else  % say == 2
    returnval.predmat = predmat_log10;
    returnval.log10A = linear_coeffs(1);
    returnval.B = -linear_coeffs(2);
end

end

