function [ intensities, xinds, yinds, combodata ] = convMax( these_maxima, conv_map )
% If no hexagonal window masks, supply BW from imregionalmax as
% these_maxima.
%
% Nathanael Kazmierczak

lininds_maxima = find(these_maxima);
[xinds,yinds] = ind2sub(size(these_maxima),lininds_maxima);
intensities = conv_map(lininds_maxima);
% Note that combo data doesn't necessarily have to be 
combodata = [intensities,xinds,yinds];

end

