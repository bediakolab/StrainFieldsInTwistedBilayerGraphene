function [ A_trimmed, rowindices, colindices ] = trimArray( A,trim_width,B )
% Nathanael Kazmierczak, 03/2020
%
% B is an array to replace the trim portion with.

if nargin < 3
    B = [];
end

[r,c] = size(A);
rlower = trim_width+1;
rupper = r - trim_width;
clower = trim_width+1;
cupper = c - trim_width;

if isempty(B)
    if r == 1
        A_trimmed = A(1,clower:cupper);
        rowindices = 1;
        colindices = clower:cupper;
    elseif c == 1
        A_trimmed = A(rlower:rupper,1);
        rowindices = rlower:rupper;
        colindices = 1;
    else
        A_trimmed = A(rlower:rupper,clower:cupper);
        rowindices = rlower:rupper;
        colindices = clower:cupper;
    end
else
    A_trimmed = A;
    A_trimmed(1:rlower,:) = B(1:rlower,:);
    A_trimmed(rupper:end,:) = B(rupper:end,:);
    A_trimmed(:,1:clower) = B(:,1:clower);
    A_trimmed(:,cupper:end) = B(:,cupper:end);
end


end

