function [strain_struct] = computeStrain(unstrained_lattice,lattice_storage)
% The entire nD array of lattice vectors should be passed here (for one
% graphene at a time). They will be packaged up into a spatial map of the
% same size as the third.
%
% Return as a strain struct so that there are not eight variables kicking
% around in our workspace.

storsize = size(lattice_storage);
if length(storsize) < 4
    storsize(4) = 1;
end

strain_struct.exx_map = zeros(storsize(3),storsize(4));
strain_struct.eyy_map = zeros(storsize(3),storsize(4));
strain_struct.exy_map = zeros(storsize(3),storsize(4));
strain_struct.theta_map = zeros(storsize(3),storsize(4));

for i = 1:storsize(3)
    for j = 1:storsize(4)
        this_strained_lattice = lattice_storage(:,:,i,j);
        if any(isnan(this_strained_lattice))
            strain_struct.exx_map(i,j) = nan;
            strain_struct.eyy_map(i,j) = nan;
            strain_struct.exy_map(i,j) = nan;
            strain_struct.theta_map(i,j) = nan;
            continue
        end
        this_strained_lattice = this_strained_lattice(2:3,:);  % getting rid of center coords which we care about not at all here.
        T = (unstrained_lattice'\this_strained_lattice')';
        strain_struct.exx_map(i,j) = T(1,1) - 1;
        strain_struct.eyy_map(i,j) = T(2,2) - 1;
        strain_struct.exy_map(i,j) = T(1,2) + T(2,1);
        strain_struct.theta_map(i,j) = T(2,1) - T(1,2);
    end
end



end

