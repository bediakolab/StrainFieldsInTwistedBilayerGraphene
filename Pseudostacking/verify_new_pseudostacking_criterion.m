% verify_new_pseudostacking_criterion.m
%
% To test and visualize the modified Wiger-Seitz partition of the
% displacement half-hexagon.
%
% This script also gives the limiting values expected in the case of no
% reconstruction. Empirically, it seems we expect about 40% AA, 30% AB, and
% 30% SP under this criterion, but would have to work this out analytically
% to be sure.
%
% Nathanael Kazmierczak, 07/11/2020

ub = 3;
lb = -3;
N = 300;
xbase = linspace(lb,ub,N);
ybase = linspace(lb,ub,N);
[xdisps,ydisps] = meshgrid(xbase,ybase);
dispfield = cat(3,xdisps,ydisps);

% Set this up to run on the incoming displacement raster generated
% manually.
make_plots = 1;
r = 'modified wigner-seitz';
filterstruct = [];
use_annealed = 3;  % The toggle for the manual input
croprange = -1;
m4 = FourDSTEM_Analysis_Engine;
[AA_perc,AB_perc,SP_perc,pseudostacking_mat] = m4.getPseudostackingAssignments(make_plots,r,filterstruct,use_annealed,croprange,dispfield);
figure; 
imagesc([lb,ub],[lb,ub],pseudostacking_mat);
axis equal
set(gca,'ydir','normal');
colormap(parula);
colorbar;
hold on;
plotFullDisplacementHexagons(gca);
