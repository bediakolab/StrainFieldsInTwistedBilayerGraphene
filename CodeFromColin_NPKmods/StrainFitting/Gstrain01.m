function [sFit] = Gstrain01(stack4D,CBEDvacuum)

sFit.stackSize = size(stack4D);
sFit.CBEDmean = mean(mean(stack4D,3),4);
sFit.CBEDmax = max(max(stack4D,[],3),[],4);
sFit.CBEDvacuum = CBEDvacuum;

end