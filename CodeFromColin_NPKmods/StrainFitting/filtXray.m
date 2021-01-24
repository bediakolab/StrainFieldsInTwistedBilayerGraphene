function [imageOut] = filtXray(imageIn,thresh)

imageMed = medfilt2(imageIn,[3 3],'symmetric');
sub = abs(imageIn - imageMed) > thresh;
imageOut = imageIn;
imageOut(sub) = imageMed(sub);
end