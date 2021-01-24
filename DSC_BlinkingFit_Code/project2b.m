function [ magnitudes ] = project2b( vec,mat )
%Projects a vector onto each of the vectors in a matrix.
% Assume they are all row vectors.

% v = repmat(vec,size(mat,1),1)./repmat(norm(vec)^2,size(mat,1),2);

% v = mat./repmat(sum(mat.^2,2),1,2);  % because we need the squared magnitude of each vector
% vec*mat'; % the inner products
% p = repmat(v,1,size(vec,1)).*repmat((mat*vec'),1,2);  % inner products can remain the same, they commute
magnitudes = mat*vec';

end

