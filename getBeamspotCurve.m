function [frac_overlap_storage] = getBeamspotCurve(beam_radius,reconstruction_size)
% Settings for the lattice are always going to be constant here

RESOLUTION = 0.02;

xbase = -2:RESOLUTION:10;
ybase = -2:RESOLUTION:8;
[xspace,yspace] = meshgrid(xbase,ybase);

xvalues = 2:RESOLUTION:8;  % For the beamspot
yvalue = 2*sqrt(3);


% make hex mask
x0_init = 0;
y0_init = 0;
d = 1;
xnum = 10;
ynum = 8;
[points, centers] = getHexagonalLattice(x0_init,y0_init,d,xnum,ynum);

[hexmask] = generateHexagonalLatticeMask(xspace,yspace,centers(:,1),centers(:,2),reconstruction_size);
frac_overlap_storage = zeros(numel(xvalues),1);

for i = 1:numel(xvalues)
    circlemask = isInCircle(xspace,yspace,xvalues(i),yvalue,beam_radius);
    bothmask = (circlemask & hexmask);
    totalpixels_incircle = nnz(circlemask);
    hexpixels_incircle = nnz(bothmask);
    frac_overlap_storage(i) = hexpixels_incircle/totalpixels_incircle;
end


end

