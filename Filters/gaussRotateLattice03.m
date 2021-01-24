function [ rotated_xy_lat ] = gaussRotateLattice03( ang, d, xy_lat, cellDimXY )
% Nathanael Kazmierczak, 12/13/2019
%
% Wrapper function for gaussRotateOnAA.m, translating the lattice for
% rotation around all five of the AA-stacked regions.
%
% Modified 12/19/2019 from gaussRotateLattice.m to implement an additional
% four rotations outside of the unit cell on the AA domain extending out 
% from each edge. This will make the reconstructed Moiré cell have flat 
% sides, which is not currently the case.
%
% Modified 04/04/2020 to add eight additional rotations to make the
% AA regions the same size when doing modulo arithmetic.

x_shift = 0;
y_shift = 0;
[ xy_lat ] = gaussRotateOnAA( ang, d, xy_lat, x_shift, y_shift );
x_shift = cellDimXY(1);
y_shift = 0;
[ xy_lat ] = gaussRotateOnAA( ang, d, xy_lat, x_shift, y_shift );
x_shift = 0;
y_shift = cellDimXY(2);
[ xy_lat ] = gaussRotateOnAA( ang, d, xy_lat, x_shift, y_shift );
x_shift = cellDimXY(1);
y_shift = cellDimXY(2);
[ xy_lat ] = gaussRotateOnAA( ang, d, xy_lat, x_shift, y_shift );
x_shift = cellDimXY(1)/2;
y_shift = cellDimXY(2)/2;
[ xy_lat ] = gaussRotateOnAA( ang, d, xy_lat, x_shift, y_shift );

% The four extra rotations
x_shift = cellDimXY(1)/2;
y_shift = -cellDimXY(2)/2;
[ xy_lat ] = gaussRotateOnAA( ang, d, xy_lat, x_shift, y_shift );
x_shift = cellDimXY(1)/2;
y_shift = cellDimXY(2)*3/2;
[ xy_lat ] = gaussRotateOnAA( ang, d, xy_lat, x_shift, y_shift );
x_shift = -cellDimXY(1)/2;
y_shift = cellDimXY(2)/2;
[ xy_lat ] = gaussRotateOnAA( ang, d, xy_lat, x_shift, y_shift );
x_shift = cellDimXY(1)*3/2;
y_shift = cellDimXY(2)*1/2;
[ xy_lat ] = gaussRotateOnAA( ang, d, xy_lat, x_shift, y_shift );

% % The eight extra rotations
x_shift = -cellDimXY(1)/2;
y_shift = -cellDimXY(2)/2;
[ xy_lat ] = gaussRotateOnAA( ang, d, xy_lat, x_shift, y_shift );
x_shift = 0;
y_shift = -cellDimXY(2);
[ xy_lat ] = gaussRotateOnAA( ang, d, xy_lat, x_shift, y_shift );
x_shift = -cellDimXY(1)/2;
y_shift = 3*cellDimXY(2)/2;
[ xy_lat ] = gaussRotateOnAA( ang, d, xy_lat, x_shift, y_shift );
x_shift = 0;
y_shift = 2*cellDimXY(2);
[ xy_lat ] = gaussRotateOnAA( ang, d, xy_lat, x_shift, y_shift );
x_shift = cellDimXY(1);
y_shift = -cellDimXY(2);
[ xy_lat ] = gaussRotateOnAA( ang, d, xy_lat, x_shift, y_shift );
x_shift = 3*cellDimXY(1)/2;
y_shift = -cellDimXY(2)/2;
[ xy_lat ] = gaussRotateOnAA( ang, d, xy_lat, x_shift, y_shift );
x_shift = 3*cellDimXY(1)/2;
y_shift = 3*cellDimXY(2)/2;
[ xy_lat ] = gaussRotateOnAA( ang, d, xy_lat, x_shift, y_shift );
x_shift = cellDimXY(1);
y_shift = 2*cellDimXY(2);
[ xy_lat ] = gaussRotateOnAA( ang, d, xy_lat, x_shift, y_shift );




rotated_xy_lat = xy_lat;

end

