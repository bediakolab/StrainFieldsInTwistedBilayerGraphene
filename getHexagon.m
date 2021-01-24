function [points] = getHexagon(x0,y0,d)
% Finds the vertices of a regular hexagon centered at (x0,y0) and having
% length d to each edge.
points = [x0-d, y0+d/sqrt(3);
          x0-d, y0-d/sqrt(3);
          x0+d, y0+d/sqrt(3);
          x0+d, y0-d/sqrt(3);
          x0, y0+2*d/sqrt(3);
          x0, y0-2*d/sqrt(3)];
end

