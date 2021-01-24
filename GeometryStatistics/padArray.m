function [ newA ] = padArray( A, padval )
% function for duplicating the rows and columns outwards so that filtering
% algorithms don't unjustly chop off stuff near the boundary.

[r,c] = size(A);
newA = zeros(r+2*padval,c+2*padval);
newA(padval+1:end-padval,padval+1:end-padval) = A;
% Add to top horizontal face of matrix
newA(1:padval,padval+1:end-padval) = repmat(A(1,:),padval,1);
newA(1:padval,padval+1:end-padval) = repmat(A(1,:),padval,1);

end

