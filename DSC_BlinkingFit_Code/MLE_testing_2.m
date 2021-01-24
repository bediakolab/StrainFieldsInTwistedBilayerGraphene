% MLE_testing_2.m

AA_probs = ones(1,12);
AB_probs = horzcat(0.25*ones(1,6),ones(1,6));
% loglik_point = @(lam,k) -(k.*log(lam) - lam + log(factorial(k)));

niter = 10000;
ssize = 5;
for i = 1:niter
    kAA = poissrnd(repmat(AA_probs,ssize,1));
    loglik_point_fit = @(lam,filter) loglik_point(ones(size(kAA)).*lam,kAA);
end
