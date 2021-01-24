function [ reduced_zone_disps ] = extendedZoneDisp2ReducedZoneDisp( extended_zone_disps )
% For use in determining the moire angle, strain mapping against mean, etc.
%
% Assumes individual disps are row vectors. May have may row vectors
%
% Nathanael Kazmierczak, 05/05/2020


perm_val = 10;
basis = permn(-perm_val:perm_val,2);
[ v1, v2, hexagon_lattice_constant ] = getDSCBasisVectors();
w1 = v1 + v2;
w2 = 2*v2 - v1;
W = [w1',w2'];
genpoints = W*basis';

genpointsx = genpoints(1,:);
genpointsy = genpoints(2,:);
extended_zone_dispsx = extended_zone_disps(:,1);
extended_zone_dispsy = extended_zone_disps(:,2);

xdiffs = genpointsx - extended_zone_dispsx;
ydiffs = genpointsy - extended_zone_dispsy;
dists = sqrt(xdiffs.^2 + ydiffs.^2);
[~,mindistidxs] = min(dists,[],2);

w_points = genpoints(:,mindistidxs)';  % which correspond row-by-row to the displacement vector
w_translated_disps = extended_zone_disps - w_points;
reduced_zone_disps = w_translated_disps;
reduced_zone_disps(reduced_zone_disps(:,2) < 0,:) = -reduced_zone_disps(reduced_zone_disps(:,2) < 0,:);

sametest = extended_zone_disps == reduced_zone_disps;
n_unchanged = nnz(sametest(:,1) & sametest(:,2));
ntot = size(extended_zone_disps,1);
n_changed = ntot - n_unchanged;
fprintf('extendedZoneDisp2ReducedZoneDisp() modified %d of %d vectors in the displacement field.\n',n_changed,ntot);

end

