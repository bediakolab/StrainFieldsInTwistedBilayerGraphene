function [ residuals, J ] = HamishTrigFittingFunctionResidualsWrapper( variables, fitting_values, weight_vector )
% Residuals function for interferometry fitting of the 4DSTEM data. Uses
% Hamish's approach of not only fitting the displacement vectors, but also
% the magnitudes of the cosine prefactors. 
%
% To be anonymized in the script. Pass the fitting values in linearly to
% match as well. Needs to be disks 1:12 for displacement vector 1, then
% 1:12 for displacement vector 2, etc, where the displacement vectors are
% linearly opened. 
%
% Nathanael Kazmierczak, 04/06/2020

prefactors = variables(1:12)';
k = numel(variables(13:end));
displacement_guesses = reshape(variables(13:end),[2,k/2]);
if nargout > 1
    [ pred_vals, J ] = HamishTrigFittingFunction( displacement_guesses, prefactors );
%     J = -J;
else
    [ pred_vals ] = HamishTrigFittingFunction( displacement_guesses, prefactors );
end
if any(~(weight_vector == 1))
    n = numel(fitting_values);
    onesmat = ones(12,n/12).*weight_vector';
    residuals = (fitting_values - pred_vals).*onesmat(:);
    if nargout > 1
        % this would work for a proper weight vector, but wants to create a
        % full array. Doubtless there is a better way to do this later.
%         J = J.*onesmat(:);
%       In fact we can probably leave this out all together, since the
%       residuals for the skipped disks are always zero anyways. But to be
%       safe, leave it in for now.
        J(~logical(onesmat(:)),:) = 0;
    end
else
    residuals = (fitting_values - pred_vals);
end


end

