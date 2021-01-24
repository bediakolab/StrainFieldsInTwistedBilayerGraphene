function [ colormap_array ] = makeColorMapFromImage( image_filepath, top_trim, bottom_trim )
% Function for taking a screenshot of a colorbar (as in, say, Pablo or
% Philip's papers) and turning it into a colormap for use in Matlab.
%
% Nathanael Kazmierczak, 07/15/2020

if nargin < 2
    top_trim = 0;
    bottom_trim = 0;
end

[im] = imread(image_filepath);
im = im((bottom_trim+1):(end-top_trim),:,:);

R = im(:,:,1);
G = im(:,:,2);
B = im(:,:,3);
[r,c] = size(R);

Rmean = mean(R,2);
Gmean = mean(G,2);
Bmean = mean(B,2);

colormap_array = horzcat(Rmean,Gmean,Bmean)./256;

end

