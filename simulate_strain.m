% simulate_strain.m
%
% Script for probing the effects of strain on a hexagonal lattice, and in
% particular whether we can use the center of mass of the bragg disk
% positions to "hexagonalate" the position of the beam under the beamstop.

% Needs to be centered for the rotation matrix to work.
xbase = -256:256;
ybase = -256:256;
[xspace,yspace] = meshgrid(xbase,ybase);
x0 = 0;
y0 = 0;
d = 80;
[ref_lattice] = getHexagon(x0,y0,d)';
ref_lattice = orderHexagonVertices(ref_lattice')';
ref_lattice_plot = [ref_lattice, ref_lattice(:,1)];
% So that the graph is pretty with a connected hexagon
figure; scatter(ref_lattice(1,:)',ref_lattice(2,:)','filled','b');
hold on;
xlim([xbase(1),xbase(end)]);
ylim([ybase(1),ybase(end)]);

e_xx = 0.4;
e_yy = 0.3;
e_xy = 0.7;
theta = 0.9;

[strain_transform] = buildStrainTransform(e_xx,e_yy,e_xy,theta);
% Here I believe it is the case that the Ophus paper has a typo. It seems
% that they have a right-transformation, but from the jupyter notebook
% tutorial, we should have a transformation matrix acting from the left.
strained_lattice = strain_transform*ref_lattice;
strained_lattice = orderHexagonVertices(strained_lattice')';

strained_lattice = strained_lattice;

strained_lattice_plot = [strained_lattice, strained_lattice(:,1)];
hold on;
scatter(strained_lattice(1,:),strained_lattice(2,:),'r','filled');
plot(ref_lattice_plot(1,:)',ref_lattice_plot(2,:)','b-');
plot(strained_lattice_plot(1,:)',strained_lattice_plot(2,:)','r-');
% figure out where the quiver arrows should start; at the COM.
com_unstrained = mean(ref_lattice,2)
com_strained = mean(strained_lattice,2)
ax = gca;
% note that for whatever reason, the quiver plot seems to quit working when
% we displace the origin of the hexagon. But given that plot and scatter
% are fine, I think this is fine.
quiver(ax,repmat(com_unstrained(1),6,1),repmat(com_unstrained(2),6,1),ref_lattice(1,:)',ref_lattice(2,:)','b');
quiver(ax,repmat(com_strained(1),6,1),repmat(com_strained(2),6,1),strained_lattice(1,:)',strained_lattice(2,:)','r');
% Apparently if you want to restrict the legend to only two things, you
% need to have it come after all the plotting calls.
legend('Unstrained','Strained');
mystring = strcat(['Strain simulation: e_x_x=',num2str(e_xx),', e_y_y=',num2str(e_yy),', e_x_y=',num2str(e_xy),', \theta=',num2str(theta)]);
title(mystring);
xlabel('x reciprocal');
ylabel('y reciprocal');




