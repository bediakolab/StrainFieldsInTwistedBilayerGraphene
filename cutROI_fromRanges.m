function [ roicut ] = cutROI_fromRanges( DP, IRangeCut, JRangeCut )
% Simple function to use a present GUI I and J range to also make the same
% cut out of several other images, for, say, disk registration in an ROI
% region defined only once.

roicut = DP(IRangeCut,JRangeCut);

end

