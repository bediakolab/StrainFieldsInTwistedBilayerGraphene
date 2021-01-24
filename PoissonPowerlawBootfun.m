function [origin_coords] = PoissonPowerlawBootfun(lineardata)%,c_initial_MLE_guesses)
% Function for taking in data and spitting out location estimates, to be
% passed as bootfun to Matlab's bootstrp in order to ascertain the variance
% of the centerpoint estimation.
%
% Note that the quick-and-dirty way of doing this is not very principled as
% far as a bootstrap goes, because I will be resampling from masked regions
% as well and thus the sample size will not be the same for each resampling
% iteration. But nonetheless this should give a quick idea.
%
% On second thought, because bootfun resamples from the rows and not
% pointwise, I'm not going to be able to do that. Need to perform all
% calculations linearly. Cannot easily reshape into a residuals matrix.

y0_init = (max(lineardata(:,2)) - min(lineardata(:,2)))/2;
x0_init = (max(lineardata(:,3)) - min(lineardata(:,3)))/2;


log10A_init_MLE = 5;
B_init_MLE = -2;
c_initial_MLE_guesses = [x0_init,y0_init,log10A_init_MLE,B_init_MLE];
powerLawMLEPoissFit = @(c) getFullRadialPoissonMaxLikFit_bootstrap(c,lineardata);
options = optimset;
options.Display = 'iter';
options.TolFun = 1e-7;
options.TolX = 1e-7;

c_MLE_optimal = fminsearch(powerLawMLEPoissFit,c_initial_MLE_guesses,options);
origin_coords = c_MLE_optimal(1:2);

end

