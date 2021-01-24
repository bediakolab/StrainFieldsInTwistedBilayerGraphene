function [ Aout ] = bwneighborerode( A, inclusive_thresh, n, connectivity_guard_flag )
% Nathanael Kazmierczak, 04/01/2020
% Takes a binary mask
% if connectivity_guard_flag == 1, then only have separated filters
% connectivity_guard_flag == 2, then allow adjacent filters
if isinf(n)
    n = 100000;
end
A = uint8(A);

% disp('initial nnz A');
% nnz(A)
for i = 1:n
    init_nnz = nnz(A);
    neighbor_num = imfilter(uint8(A),ones(3),'conv');
    to_flip = (neighbor_num <= inclusive_thresh) & logical(A);
    if connectivity_guard_flag > 0
        new = 1;
        if new
            [r,c] = size(A);
            [ filterstor ] = getFirstOrderConnectivityFilters(connectivity_guard_flag-1);
            numfilt = numel(filterstor);
            convres_stor = false(r,c,numfilt);
            for k = 1:numfilt
                convres_stor(:,:,k) = filter2(filterstor{k},A) == 2;
            end
            linkpoints = any(convres_stor,3) & A;
            neighbors = filter2(ones(3),linkpoints);
            A(to_flip & ~linkpoints & ~neighbors) = 0;
        else
            to_flip = find(to_flip);
            i
            to_delete = [];
            conres_init = bwconncomp(A);
            init_conn_comp = conres_init.NumObjects;
            for j = 1:numel(to_flip)
                if ~mod(j,1000)
                    fprintf('%d of %d\n',j,numel(to_flip));
                end
                trialA = A;
                trialA(to_flip(j)) = 0;
                conres_final = bwconncomp(trialA);
                final_conn_comp = conres_final.NumObjects;
                if init_conn_comp == final_conn_comp
                    A = trialA;
                end  % else continue, leaving that pixel as it was.
            end
        end
    else
        A(to_flip) = 0;
    end
    new_nnz = nnz(A);
    if init_nnz == new_nnz
        disp('Neighbor filtering converged before reaching maximum iterations.')
        break
    end
end
Aout = logical(A);
% disp('final nnz A');
% nnz(Aout)

end

