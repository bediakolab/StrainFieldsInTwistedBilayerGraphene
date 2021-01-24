function [ miss_count ] = AAmask_linearfitfun( x0,y0,r,xspacecropped_org,yspacecropped_org,spacedata_lin,datavals_lin )
% Here spacedata_lin are giving the coordinates of the
% data points the fit is actually compared against.
% For bootstrapping, this will not form a grid. Lookup correspondence on a
% pixel by pixel basis using linear indexing.

[ linearmask ] = AAmask_linearpredfun( x0,y0,r,xspacecropped_org,yspacecropped_org );
% use linear indexing to pull out the predicted values corresponding to the
% rows of the datavals_lin.
predictedvals = linearmask(spacedata_lin);
misses = predictedvals-datavals_lin;
miss_count = nnz(misses);

end

