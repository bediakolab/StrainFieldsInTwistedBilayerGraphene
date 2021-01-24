function [ pred_vals ] = blinkingFittingFunction( DSC_guess )
% The pred vals give predicted blinking pattern. Anonymize the fit in
% driver.
%
% Nathanael Kazmierczak, Dec 2019.

% load('12202019FittingFunctions');  % loads func_cell containing the anonymized fitting functions
persistent func_cell_to_use

if isempty(func_cell_to_use)  % ~exist('func_cell_to_use','var') ||
    load('12232019FittingFunctionsZeroed.mat');  % populates func_cell
    func_cell_to_use = func_cell;
    %clearvars -except DSC_guess func_cell_to_use;
end

pred_vals = zeros(size(DSC_guess,1),12);
for i = 1:12
    this_func = func_cell_to_use{i};
    pred_vals(:,i) = this_func(DSC_guess);
end

end

