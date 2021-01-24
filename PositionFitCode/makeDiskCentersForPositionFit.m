function [ disk_centers_G1, disk_centers_G2, indices_G1, indices_G2, lattice_vectors_G1, lattice_vectors_G2, beam_center_coords ] = makeDiskCentersForPositionFit( diffraction_patterns )%, beam_center_coords )
% Nathanael Kazmierczak, Feb 2020
%
% Modification of makeDiskCenters to get the initial lattice for the disk
% position fitting.
%
% Assume diffraction_patterns is a 4D array such as would come out of the
% original.
%
% The definition of the disks is such where the topmost disk in the inner
% ring should be done first and proceed clockwise, while the outer ring
% should go the same way except the first outer disk should be the one
% trailing the first inner one on the clockwise rotation.

plotDP(diffraction_patterns(:,:,1,1));
fighandle = gcf;
% xlim([2.263892128279883e+02 3.706107871720116e+02]);
% ylim([2.272580174927113e+02 3.680043731778425e+02]);
totalsize = size(diffraction_patterns);
disk_centers_G1 = zeros(0,2);
disk_centers_G2 = zeros(0,2);
trial_basis = zeros(2,2);

i = 0;
while true
    i = i+1;
    while true
        fprintf('Zoom graph to prepare to click the %dth disks for lattices 1 and 2.\n',i);
        thisans = input('If the desired disk is not visible, enter 0; otherwise any entrance will continue.');
        if thisans == 0
            close(fighandle);
            disp('Displaying a randomly different diffraction pattern.');
            randind1 = randi(totalsize(3));
            randind2 = randi(totalsize(4));
            fighandle = plotDP(diffraction_patterns(:,:,randind1,randind2));
        else
            break
        end
    end
    
    
    fprintf('Click once to set the center of the %dth disk, lattice 1.\n',i);
    disk_centers_G1(end+1,:) = ginput(1);
    fprintf('Click once to set the center of the %dth disk, lattice 2.\n',i);
    disk_centers_G2(end+1,:) = ginput(1);
        
    tf = input('1/0: Are you done selecting disks for position fitting?');
    if tf
        break
    else
        figure(fighandle);
    end
    
end

fighandle = plotDP(diffraction_patterns(:,:,1,1));
hold on
% viscircles(disk_centers_G1,radius*ones(size(disk_centers_G1,1),1),'Color','g');
% viscircles(disk_centers_G2,radius*ones(size(disk_centers_G2,1),1),'Color','b');
scatter(disk_centers_G1(:,1),disk_centers_G1(:,2),'Filled','g');
scatter(disk_centers_G2(:,1),disk_centers_G2(:,2),'Filled','b');

input('Ensure graph is ready for defining trial basis and pattern center. Press any key to continue.');
fprintf('Click once to set trial basis vector 1.\n');
trial_basis(1,:) = ginput(1);
fprintf('Click once to set trial basis vector 2.\n');
trial_basis(2,:) = ginput(1);
fprintf('Click once to set the center of the diffraction spots.\n');
beam_center_coords = ginput(1);

%% Fit a hexagonal lattice to the selected points
% The following code comes from makeHexagonalLatticeMask.m
% First, index using the first point as the center and the second and third
% points as the tips of the initial guess lattice vectors.

xy0 = beam_center_coords;  % should be a column vector.
if size(xy0,2) > size(xy0,1)
    xy0 = xy0';
end
v_inits = trial_basis' - repmat(xy0,1,size(trial_basis',2));
disk_centers_G1 = disk_centers_G1';
disk_centers_G2 = disk_centers_G2';
allpoints_centered_G1 = disk_centers_G1 - repmat(xy0,1,size(disk_centers_G1,2));
allpoints_centered_G2 = disk_centers_G2 - repmat(xy0,1,size(disk_centers_G2,2));
indices_G1 = v_inits\allpoints_centered_G1;
indices_G1 = round(indices_G1);
indices_G2 = v_inits\allpoints_centered_G2;
indices_G2 = round(indices_G2);

% Next, having indexed all points, figure out which lattice vectors best
% fit those points.
lattice_vectors_G1 = allpoints_centered_G1/indices_G1;
lattice_vectors_G2 = allpoints_centered_G2/indices_G2;


end

