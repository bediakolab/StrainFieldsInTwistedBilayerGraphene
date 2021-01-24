function [ linearmask ] = AAmask_linearpredfun( x0,y0,r,xspace,yspace )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

mask = isInCircle(xspace,yspace,x0,y0,r);
linearmask = mask(:);


end

