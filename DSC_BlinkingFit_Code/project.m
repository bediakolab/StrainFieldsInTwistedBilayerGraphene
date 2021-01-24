function [ p, magnitudes ] = project( mat,vec )
%Projects matrix mat onto vec.
% Assume they are all row vectors.

v = repmat(vec,size(mat,1),1)./repmat(norm(vec)^2,size(mat,1),2);
% vec*mat'; % the inner products
p = v.*repmat((mat*vec'),1,2);
magnitudes = mat*vec';

end

