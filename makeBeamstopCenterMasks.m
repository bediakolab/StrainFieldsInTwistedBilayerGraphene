function [beamstop_masks,Irng,Jrng,this_maskeddata] = makeBeamstopCenterMasks(data,use_teeny_mask,teeny_mask_radius_factor,draw_hard_delete_region,THRESHOLD)
% Function for masking beamstop in preparation for a power law fit.
% Wrapper for code originally developed in
% radial_powerlaw_fit_beamcenter_determination.m
% Nathanael Kazmierczak, 08/11/2019
% Bediako Lab, UC Berkeley
%
% Outputs:
% maskeddata is what powerlaws should be fit to (e.g. Poisson MLE).
% Irng, Jrng give the coordinates of the extracted region in first matrix
% dimension (y) and second matrix dimension (x) respectively. Needed for
% conversion to global coordinates later.
%
% While the idea here is to generate masks that can then be applied on
% later DPs without the need to draw, the process of ascertaining the
% mask's validity is bound up with application to a particular DP, so that
% masked data will be returned in case the user wants it.


%% Get the inner circle region

%THRESHOLD = 50;
% THRESHOLD = 100;
if nargin < 6
    THRESHOLD = 100;
end

% figure; image(uint8(data),'CDataMapping','scaled');
plotDP(data);
if ~verLessThan('matlab','9.2')
    disp('Please draw a mask to demarcate the center of the image for power law fitting.');
    mycircle = drawcircle;
    disp('Finished drawing circular mask for power law fit.');
    circle_mask = createMask(mycircle);
else
    disp('Please click twice on the graph, once to position the center of the circle, the second time to position the radius.');
    [xc,yc] = ginput(2);
    r = sqrt((xc(1)-xc(2))^2 + (yc(1)-yc(2))^2);
    xbase = 1:size(data,1);
    ybase = 1:size(data,2);
    [xspace,yspace] = meshgrid(xbase,ybase);
    circle_mask = isInCircle(xspace,yspace,xc(1),yc(2),r);
end

% Now make a new plot.
maskeddata = double(data);
maskeddata(~circle_mask) = -1;
[Icoords,Jcoords] = ind2sub(size(circle_mask),find(circle_mask));
Irng = [min(Icoords),max(Icoords)];
Jrng = [min(Jcoords),max(Jcoords)];
maskeddata = maskeddata(Irng(1):Irng(2),Jrng(1):Jrng(2));
newdims = size(maskeddata);
org_radius = newdims(1)/2;
org_circle_mask_data = maskeddata;
figure;
image(uint8(maskeddata),'CDataMapping','scaled');

%% Indicate the region of the beamstop
% Anything below a threshold in this region will be set to nan, so try to
% avoid hitting off-target regions.
disp('Please draw a polygon enclosing the region of the beamstop but not filtering out regions for the power law fit.')
if ~verLessThan('matlab','9.2')
    mypolygon = drawpolygon;
    beamstop_mask = createMask(mypolygon);
else
    beamstop_mask = roipoly;
end
maskeddata(beamstop_mask & maskeddata < THRESHOLD) = -1;  % let this be the sign to ignore
% The first delete above is a 'soft' delete because of the threshold.
% To prevent rollover without etching too many outer pixels, futhermore 
% have a 'hard' delete without a threshold for half of the circle radius.
% Define an inner circular mask by halving the radius on the initial circle
% roi.
% Expand the beamstop mask by three pixels in each direction to ensure we
% are not getting data that is rolling over. Employ an iterative boundary
% building using Matlab's nice functions.
bmask = (maskeddata == -1);
bmask2 = boundarymask(bmask);
bmask3 = boundarymask(bmask2);
%bmask4 = boundarymask(bmask3);
beamstop_boundary_mask = bmask | bmask2 | bmask3;

% in what region will we apply hard delete of the thickened boundary?
if draw_hard_delete_region == 1
    % In this case it is a polygon and not a circle, lol.
    disp('Please draw a polygon specifying the inner hard delete region.');
    mycircle2 = drawpolygon;
    inner_circle_mask = createMask(mycircle2);
elseif draw_hard_delete_region == 2
    inner_circle_mask = beamstop_mask;  % hard delete everything (02/02/2020)
    maskeddata(inner_circle_mask) = -1;
else
    mycircle2 = mycircle;
    mycircle2.Radius = org_radius / 2;
    inner_circle_mask = createMask(mycircle2);
    inner_circle_mask = inner_circle_mask(Irng(1):Irng(2),Jrng(1):Jrng(2));
end
maskeddata(inner_circle_mask & beamstop_boundary_mask) = -1;

if use_teeny_mask
    % Create also a smallest inner circle mask for deletion, citing concerns
    % about elliptical distortions.
    mycircle3 = mycircle;
    mycircle3.Radius = org_radius*teeny_mask_radius_factor;
    teeny_circle_mask = createMask(mycircle3);
    teeny_circle_mask = teeny_circle_mask(Irng(1):Irng(2),Jrng(1):Jrng(2));
    maskeddata(teeny_circle_mask) = -1;  % for exclusion from the fit.
end

figure;
image(uint8(maskeddata),'CDataMapping','scaled');
title 'Masked data from makeBeamstopCenterMasks.m'
% Now, this data is good enough to do the actual power law fits.

% Return data as a struct, because teeny mask may not be present and there
% a lot of masks, so a bit cumbersome.
beamstop_masks.circle_mask = circle_mask;
beamstop_masks.beamstop_mask = beamstop_mask;
beamstop_masks.inner_beamstop_mask = beamstop_boundary_mask;
beamstop_masks.inner_circle_mask = inner_circle_mask;
if use_teeny_mask
    beamstop_masks.teeny_circle_mask = teeny_circle_mask;
end
beamstop_masks.Irng = Irng;
beamstop_masks.Jrng = Jrng;
this_maskeddata = maskeddata;

end

