function [ residuals ] = trigFittingFunction2outer( params, disk_average_storage )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

prefactors = params(1:12)';
DSC_params = params(13:end);
DSC_params_tofit = reshape(DSC_params,[2,numel(DSC_params)/2])';

% assemble DSC and return result. 
[s1,s2,~] = size(disk_average_storage);
[ pred_vals ] = trigFittingFunction2( prefactors, DSC_params_tofit );
% pred_vals(isnan(pred_vals)) = 0;
% disk_average_storage(isnan(disk_average_storage)) = 0;
residuals = pred_vals(:) - disk_average_storage(:);
residuals(isnan(residuals)) = 0;

end

