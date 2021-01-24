function [global_coords] = convertMaskedCoordsToOriginalCoords(masked_coords,Irng,Jrng)
% here Irng is the first dimension (i.e. y) coordinates that the masked
% matrix occupies, while Jrng (i.e. x) is the second dimension coordinates
% that the masked matrix occupies.
%
% NK 08/11/2019: This function could really use more testing to ensure I
% have the index manipulation right

% masked_coords is talking about optimization outputs, and so x is the
% first coord, to be used with Jrng
global_coords(1) = masked_coords(1) - 1 + Jrng(1);
global_coords(2) = masked_coords(2) - 1 + Irng(1);

end

