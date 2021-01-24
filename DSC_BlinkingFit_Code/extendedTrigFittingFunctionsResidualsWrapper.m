function [ residual_vals, J ] = extendedTrigFittingFunctionsResidualsWrapper( DSC_guess, scaling_constant, bound_handling_flag, third_fourth_disk_flag, fitdata, weight_vector )
% Nathanael Kazmierczak
% 03/27/2020
% Needs to be anonymized in caller.
% Applies weight vector in here so be sure not to double count.

if nargout > 1
    [ pred_vals, J ] = extendedTrigFittingFunctions( DSC_guess, scaling_constant, bound_handling_flag, third_fourth_disk_flag );
else
    [ pred_vals ] = extendedTrigFittingFunctions( DSC_guess, scaling_constant, bound_handling_flag, third_fourth_disk_flag );
end
residual_vals = fitdata - pred_vals;
residual_vals = residual_vals.*weight_vector;


end

