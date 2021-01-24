% displacement_density_plot.m
%
% For DS4S1. Try versions based on ksdensity, interpolating into a contour
% map, e.g.
%
% Nathanael Kazmierczak, 10/26/2020

if ~exist('m4','var')
addpath('/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/StrainData/AnnealedDatasetsForUpload/ChiralFixed');
load('05312020_chiralfixed_DS4_postannealing.mat');
end
load('demo_colormap.mat');

xcomps = m4.DSC_fit_storage(:,:,1);
ycomps = m4.DSC_fit_storage(:,:,2);
xcomps = xcomps(:);
ycomps = ycomps(:);
coords = [xcomps,ycomps];
figure;
inc = 0.05;
xvec = -1.4:inc:1.4;
yvec = 0:inc:1.4;
[ptsx,ptsy] = meshgrid(xvec,yvec);
ptsxn = ptsx(:);
ptsyn = ptsy(:);
fullpts = [ptsxn,ptsyn];
% ksdensity(coords,fullpts);
vals = ksdensity(coords,fullpts);
matvals = reshape(vals,size(ptsx));

% need area of half hexagon in order to know how to plot this. The scaling
% of ksdensity is automatically to unit area I believe.
[ v1, v2, hexagon_lattice_constant ] = getDSCBasisVectors();
height1 = v2(2);
height2 = v1(2) - v2(2);
len = v2(1)*2;
rectarea = height1*len;
trianglearea = 0.5*height2*len;
halfhexarea = rectarea + trianglearea;
uniformdensity = 1/halfhexarea

surf(xvec,yvec,matvals/uniformdensity);
shading interp;
% plotFullDisplacementHexagons(gca);
xlabel('x displacement (A)');
ylabel('y displacement (A)');
zlabel('Density relative to uniform distribution');
scaleColorMap( cmap, 1 );



% Alternatively, make this a flat contour plot
% contourf(xvec,yvec,matvals/uniformdensity);  % for the coarse version
contourf(xvec,yvec,matvals/uniformdensity,50,'LineStyle','none');
xlabel('x displacement (A)');
ylabel('y displacement (A)');
% zlabel('Density relative to uniform distribution');
hold on
plotFullDisplacementHexagons(gca);
scaleColorMap( cmap, 1 );
buf = 0.05;
xlim([-len/2-buf,len/2+buf]);
ylim([-buf,v1(2)+buf]);
ch = colorbar;
ylabel(ch,'Density relative to uniform distribution');

% Alternatively, 2D hist
figure;
hist3(coords);
xlabel('x displacement (A)');
ylabel('y displacement (A)');
zlabel('Density relative to uniform distribution');
