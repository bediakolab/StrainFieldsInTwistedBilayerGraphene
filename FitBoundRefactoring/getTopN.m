function [ tops ] = getTopN( array,n )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

array = fliplr(sort(array));
tops = array(1:n);


end

