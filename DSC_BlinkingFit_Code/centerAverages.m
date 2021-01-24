function [ averages_c ] = centerAverages( averages )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

% averages = averages - min(averages);
averages_c = averages./repmat(max(averages),size(averages,1),1);

end

