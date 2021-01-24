function [ path_pixels, distance ] = bwshortestpath( A, B, thinflag, inclusiveflag )
% Computes the shortest path between the pixel clumps in A and B, both
% binary images.
%
% If thinflag is true, then this will return a unique path by thinning down
% the connection matrix.
%
% A and B should each contain only a single connected component.
%
% See https://blogs.mathworks.com/steve/2011/11/01/exploring-shortest-paths-part-1/
%
% Nathanael Kazmierczak, 04/03/2020

TOL = 32;
D1 = bwdist(A, 'quasi-euclidean');
D2 = bwdist(B, 'quasi-euclidean');
D_sum = D1 + D2;
D_sum = round(D_sum*TOL)/TOL;
path_pixels = imregionalmin(D_sum);
distance = D_sum(find(path_pixels,1));
if thinflag
    path_pixels = bwmorph(path_pixels,'thin',inf);
end
if inclusiveflag
    path_pixels = path_pixels | A | B;
end
    

end

