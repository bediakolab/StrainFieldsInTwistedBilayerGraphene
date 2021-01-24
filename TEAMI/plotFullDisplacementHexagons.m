function plotFullDisplacementHexagons( axh, line_thickness )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

if nargin < 2
    line_thickness = 1;
end


[ v1, v2, hexagon_lattice_constant ] = getDSCBasisVectors();
% t = 60/180*pi;
% rotmat = [cos(t) sin(t); -sin(t) cos(t)];
% v1 = [0,hexagon_lattice_constant/sqrt(3)];
% v2 = v1*rotmat';
% vl1 = vectorLine(v2-v1,v1);
% vl2 = vectorLine(v1,v2);
% vl3 = vectorLine(v2,-v1);
% vl4 = vectorLine(v2-v1,-v1);
% vl5 = vectorLine(v1,-v2);
% vl6 = vectorLine(-v2,v1);

axes(axh);
figh = gcf;

hold on
axis equal
plotLine( figh, v1, v2, 'b', line_thickness );
plotLine( figh, v2, v2-v1, 'b', line_thickness );
plotLine( figh, -v1, v2-v1, 'b', line_thickness );
plotLine( figh, -v1, -v2, 'b', line_thickness );
plotLine( figh, v1-v2, -v2, 'b', line_thickness );
plotLine( figh, v1-v2, v1, 'b', line_thickness );

% [ pred_vals ] = trigFittingFunctions( DSC_guess, scaling_constant, bound_handling_flag );

plotLine(figh, 2*v2, 3*v2 - v1, 'b', line_thickness);
plotLine(figh, 3*v2 - 2*v1, 3*v2 - v1, 'b', line_thickness);
plotLine(figh, 3*v2 - 2*v1, 2*v2 - 2*v1, 'b', line_thickness);
plotLine(figh, 2*v2 - 2*v1, 2*v2 - 2*v1, 'b', line_thickness);
plotLine(figh, 2*v2 - 2*v1, v2-v1, 'b', line_thickness);

plotLine(figh, v1, 2*v1, 'b', line_thickness);
plotLine(figh, 2*v1 + v2, 2*v1, 'b', line_thickness);
plotLine(figh, 2*v1 + v2, 2*v2 + v1, 'b', line_thickness);
plotLine(figh, 2*v2, 2*v2 + v1, 'b', line_thickness);
plotLine(figh, 2*v2, v2, 'b', line_thickness);

plotLine(figh, v1 - 3*v1, 2*v1 - 3*v1, 'b', line_thickness);
plotLine(figh, 2*v1 + v2 - 3*v1, 2*v1 - 3*v1, 'b', line_thickness);
plotLine(figh, 2*v1 + v2 - 3*v1, 2*v2 + v1 - 3*v1, 'b', line_thickness);
plotLine(figh, 2*v2 - 3*v1, 2*v2 + v1 - 3*v1, 'b', line_thickness);
plotLine(figh, 2*v2 - 3*v1, v2 - 3*v1, 'b', line_thickness);
plotLine(figh, -2*v1, -3*v1 + v2, 'b', line_thickness);

plotLine(figh, v1 - 3*v2, 2*v1 - 3*v2, 'b', line_thickness);
plotLine(figh, 2*v1 + v2 - 3*v2, 2*v1 - 3*v2, 'b', line_thickness);
plotLine(figh, 2*v1 + v2 - 3*v2, 2*v2 + v1 - 3*v2, 'b', line_thickness);
plotLine(figh, 2*v2 - 3*v2, 2*v2 + v1 - 3*v2, 'b', line_thickness);
plotLine(figh, 2*v2 - 3*v2, v2 - 3*v2, 'b', line_thickness);
plotLine(figh, -2*v2, -3*v2 + v1, 'b', line_thickness);

plotLine(figh, v1 - 3*v2 + v1 + v2, 2*v1 - 3*v2 + v1 + v2, 'b', line_thickness);
plotLine(figh, 2*v1 + v2 - 3*v2 + v1 + v2, 2*v1 - 3*v2 + v1 + v2, 'b', line_thickness);
plotLine(figh, 2*v1 + v2 - 3*v2 + v1 + v2, 2*v2 + v1 - 3*v2 + v1 + v2, 'b', line_thickness);
plotLine(figh, 2*v2 - 3*v2 + v1 + v2, 2*v2 + v1 - 3*v2 + v1 + v2, 'b', line_thickness);
plotLine(figh, 2*v2 - 3*v2 + v1 + v2, v2 - 3*v2 + v1 + v2, 'b', line_thickness);
plotLine(figh, -2*v2 + v1 + v2, -3*v2 + v1 + v1 + v2, 'b', line_thickness);

plotLine(figh, v1 - 3*v2 + v1 + v2 - 3*v1, 2*v1 - 3*v2 + v1 + v2 - 3*v1, 'b', line_thickness);
plotLine(figh, 2*v1 + v2 - 3*v2 + v1 + v2 - 3*v1, 2*v1 - 3*v2 + v1 + v2 - 3*v1, 'b', line_thickness);
plotLine(figh, 2*v1 + v2 - 3*v2 + v1 + v2 - 3*v1, 2*v2 + v1 - 3*v2 + v1 + v2 - 3*v1, 'b', line_thickness);
plotLine(figh, 2*v2 - 3*v2 + v1 + v2 - 3*v1, 2*v2 + v1 - 3*v2 + v1 + v2 - 3*v1, 'b', line_thickness);
plotLine(figh, 2*v2 - 3*v2 + v1 + v2 - 3*v1, v2 - 3*v2 + v1 + v2 - 3*v1, 'b', line_thickness);
plotLine(figh, -2*v2 + v1 + v2 - 3*v1, -3*v2 + v1 + v1 + v2 - 3*v1, 'b', line_thickness);

plotLine( figh, v1, v2, 'r', line_thickness );
plotLine( figh, v2, v2-0.5*v1, 'r', line_thickness );
plotLine( figh, v2-0.5*v1, -v2+0.5*v1, 'r', line_thickness );
plotLine( figh, -v2+0.5*v1, -v2+v1, 'r', line_thickness );
plotLine( figh, v1-v2, v1, 'r', line_thickness );

end

