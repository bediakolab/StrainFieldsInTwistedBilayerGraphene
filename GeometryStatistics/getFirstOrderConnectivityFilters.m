function [ filtermat_stor ] = getFirstOrderConnectivityFilters( allow_adjacents_flag )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

[M] = permn(1:9, 2);
M(M(:,1)==M(:,2),:) = [];
M(M(:,1)==5 | M(:,2)==5,:) = [];

filtermat_stor = cell(0,1);
for i = 1:size(M,1)
    tp = M(i,:);
    A = zeros(3);
    A(tp) = 1;
    %     A
    if allow_adjacents_flag
        filtermat = 3*ones(3);
        filtermat(logical(A)) = 1;
        filtermat(5) = 0;
        filtermat_stor{end+1} = filtermat;
    else
        screenmat = 0.5*ones(3);
        screenmat(2,2) = 1;
        convres = filter2(screenmat,A);
        if max(convres) <= 1
            filtermat = 3*ones(3);
            filtermat(logical(A)) = 1;
            filtermat(5) = 0;
            filtermat_stor{end+1} = filtermat;
        end
    end
end


end

