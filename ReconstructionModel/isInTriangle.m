function [ mask ] = isInTriangle( xspace, yspace, vertices )
% Function for building a triangular mask, based off of the three vertices
% specified. Used in reconstructionFunctionNPK.m to give the AB rotation
% field.
%
% Vertices should be a 3x2 array, proceeding from left to right on the
% vertices.
%
% Nathanael Kazmierczak, 06/08/2020

vertices = sortrows(vertices);

v1 = vertices(1,:);
v2 = vertices(2,:);
v3 = vertices(3,:);

left_line = vectorLine(v2-v1,v1);
right_line = vectorLine(v3-v2,v2);
long_line = vectorLine(v3-v1,v1);

space = [xspace(:),yspace(:)];
if left_line > v3 
    c1 = left_line >= space; 
else
    c1 = left_line <= space;
end

if right_line > v1
    c2 = right_line >= space;
else
    c2 = right_line <= space;  % empirical
end

if long_line > v2
    c3 = long_line >= space;
else
    c3 = long_line <= space;
end

c1 = reshape(c1,size(xspace));
c2 = reshape(c2,size(xspace));
c3 = reshape(c3,size(xspace));
mask = c1 & c2 & c3;
% mask = reshape(c,size(xspace));

end

