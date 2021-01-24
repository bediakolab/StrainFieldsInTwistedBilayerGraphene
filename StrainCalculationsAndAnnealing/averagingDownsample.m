function [ newA ] = averagingDownsample( A, int_factor )
% Nathanael Kazmierczak
%
% Why on earth does Matlab not already have a built in function for this??

[r,c] = size(A);
front_edge_r = 1:int_factor:r;
back_edge_r = [front_edge_r(2:end)-1,r];
nbins_r = numel(front_edge_r);

front_edge_c = 1:int_factor:c;
back_edge_c = [front_edge_c(2:end)-1,c];
nbins_c = numel(front_edge_c);

newA = zeros(nbins_r,nbins_c);
for i = 1:nbins_r
    this_front_r = front_edge_r(i);
    this_back_r = back_edge_r(i);
    for j = 1:nbins_c
        this_front_c = front_edge_c(j);
        this_back_c = back_edge_c(j);
        newA(i,j) = mean(mean(A(this_front_r:this_back_r,this_front_c:this_back_c)));
    end
end


end

