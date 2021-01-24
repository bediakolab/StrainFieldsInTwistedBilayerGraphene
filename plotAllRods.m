function plotAllRods( DProi_storage_final,disk_storage,I,J )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

DP = DProi_storage_final{I,J};

store = disk_storage;
store(store(:,4) ~= I | store(:,5) ~= J,:) = [];
plotDPandDisks( DP, store(:,2:3) )
title(sprintf('I = %d, J = %d disk registration',I,J));

end

