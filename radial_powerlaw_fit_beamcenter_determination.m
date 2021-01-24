% radial_powerlaw_fit_beamcenter_determination.m
%
% Script attempting to find the center of the beam in the presence of a
% beamstop without using the diffraction spots.
%
% Nathanael Kazmierczak, 08/09/2019

%% Load diffraction data as normal.

addpath /Volumes/Lexar
info = h5info('3__70x70_ss4nm_0p05s_spot9_alpha=0p55_bin=4_CL=130_80kV_C2=40um__cropped.h5');
h5disp('3__70x70_ss4nm_0p05s_spot9_alpha=0p55_bin=4_CL=130_80kV_C2=40um__cropped.h5','/4DSTEM_experiment/data/datacubes/datacube_0/data');

% Attempt to read in exactly three diffraction patterns
start = [1,1,1,1];
count = [512,512,1,1];
data = h5read('3__70x70_ss4nm_0p05s_spot9_alpha=0p55_bin=4_CL=130_80kV_C2=40um__cropped.h5','/4DSTEM_experiment/data/datacubes/datacube_0/data',start,count);

% One question to answer is whether the values at the center of the beam
% are actually saturated, or whether the image display algorithms only make
% it appear to be so.
figure; surf(data); shading flat;
% Indeed, it appears that the image really is 

figure; image(uint8(data),'CDataMapping','scaled');

%% Mask the center region to be fit.
disp('Please draw a mask to demarcate the center of the image for power law fitting.');
mycircle = drawcircle;
disp('Finished drawing circular mask for power law fit.');
circle_mask = createMask(mycircle);

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
mypolygon = drawpolygon;
beamstop_mask = createMask(mypolygon);
threshold = 50;
maskeddata(beamstop_mask & maskeddata < threshold) = -1;  % let this be the sign to ignore
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
bmask4 = boundarymask(bmask3);
inner_beamstop_mask = bmask | bmask2 | bmask3;
mycircle2 = mycircle;
mycircle2.Radius = org_radius / 2;
inner_circle_mask = createMask(mycircle2);
inner_circle_mask = inner_circle_mask(Irng(1):Irng(2),Jrng(1):Jrng(2));
maskeddata(inner_circle_mask & inner_beamstop_mask) = -1;

teenymask = 0;
if teenymask
    % Create also a smallest inner circle mask for deletion, citing concerns
    % about elliptical distortions.
    mycircle3 = mycircle;
    mycircle3.Radius = org_radius/10;
    teeny_circle_mask = createMask(mycircle3);
    teeny_circle_mask = teeny_circle_mask(Irng(1):Irng(2),Jrng(1):Jrng(2));
    maskeddata(teeny_circle_mask) = -1;  % for exclusion from the fit.
end

figure;
image(uint8(maskeddata),'CDataMapping','scaled');
% Now, this data is good enough to do the actual power law fits.

%% Fit version 1: all directions, real space
% Definition of fitting functions:
% [predmat] = getFullRadialPowerLawPred(c,matsize);
% [rms_deviation] = getFullRadialPowerLawFit(c,fitdata);
powerLawFit = @(c) getFullRadialPowerLawFit(c,maskeddata);
y0_init = round(mean([1,size(maskeddata,1)]));
x0_init = round(mean([1,size(maskeddata,2)]));
log10A_init = 15;
B_init = -4;
c_initial_guesses = [x0_init,y0_init,log10A_init,B_init];
options = optimset;
options.Display = 'iter';
options.TolFun = 1e-7;
options.TolX = 1e-7;

c_optimal = fminsearch(powerLawFit,c_initial_guesses,options);

optimal_fit = getFullRadialPowerLawPred(c_optimal,size(maskeddata));
figure; subplot(1,3,1);
nan_masked_data = maskeddata;
nan_masked_data(maskeddata == -1) = nan;
surf(nan_masked_data); shading flat;
subplot(1,3,2);
surf(log10(optimal_fit)); shading flat;
subplot(1,3,3);
residuals = maskeddata - optimal_fit;
residuals(maskeddata == -1) = nan;
surf(residuals); shading flat;
colormap(flipud(bone));

% residuals_tosurf = residuals
% residuals_tosurf(isnan(residuals_tosurf)) = 0
% residuals_tosurf(residuals_tosurf <= 0) = 0.1;
% figure; surf(log10(residuals)); shading flat;
% figure; surf(log10(residuals_tosurf)); shading flat;


%% Fit version 2: all directions, log space
% This is probably the best way of doing this for two reasons: the errors
% look like they will be better distributed in log space, and the
% linearization will give us a separable linear least squares problem,
% reducing it to only two parameters (the origin coordinates) that need to
% be nonlinearly fit.

%[returnval_outer] = getFullRadialPowerLawFit_logspace(c,fitdata,optimize)
%[returnval] = getFullRadialPowerLawPred_logspace(c,matsize,fitdata,optimize)
optimize_flag = 1;
powerLawFit_logspace = @(c) getFullRadialPowerLawFit_logspace(c,maskeddata,optimize_flag);
c_init_logopt_guesses = [x0_init,y0_init];

c_logopt_optimal = fminsearch(powerLawFit_logspace,c_init_logopt_guesses,options);

optimize_flag = 0;
% remember that the fitdata should be logarithmic!!!
fitdata = maskeddata;
fitdata(fitdata <= 0) = 0.1;
[returnval] = getFullRadialPowerLawPred_logspace(c_logopt_optimal,size(maskeddata),log10(fitdata),optimize_flag);
optimal_logspace_fit = returnval.predmat;
optimal_log10A = returnval.log10A
optimal_B = returnval.B

figure; subplot(1,3,1);
nan_masked_data = maskeddata;
nan_masked_data(maskeddata == -1) = nan;
surf(nan_masked_data); shading flat;
subplot(1,3,2);
% this guy is not in logspace
surf(optimal_logspace_fit); shading flat;
subplot(1,3,3);
to_plot_resid_data = maskeddata;
to_plot_resid_data(to_plot_resid_data <= 0) = 0.1;
residuals2 = log10(to_plot_resid_data) - optimal_logspace_fit;
residuals2(maskeddata == -1) = nan;
surf(residuals2); shading flat;
colormap(flipud(copper));

%% Fit version 3: use the Poisson estimator
% The radial power law prediction mechanics still work in exactly the same
% way, so we can still use getFullRadialPowerLawPred.m in real space for
% this.
% However, now the fit function changes as it is MLE and not LSQ.
% This is fit function: [log_lik] = getFullRadialPoissonMaxLikFit(c_initial_guesses,fitdata);
log10A_init_MLE = 5;
B_init_MLE = -2;
%c_initial_MLE_guesses = [x0_init,y0_init,log10A_init_MLE,B_init_MLE];
c_initial_MLE_guesses = c_optimal;  % owing to sensitivity to initial guesses, feed the optimal from realspace opt as first guess.
powerLawMLEPoissFit = @(c) getFullRadialPoissonMaxLikFit(c,maskeddata);

c_MLE_optimal = fminsearch(powerLawMLEPoissFit,c_initial_MLE_guesses,options);

optimal_MLE_fit = getFullRadialPowerLawPred(c_MLE_optimal,size(maskeddata));
figure; subplot(1,3,1);
nan_masked_data = maskeddata;
nan_masked_data(maskeddata == -1) = nan;
surf(nan_masked_data); shading flat;
title 'Data for MLE Poiss fit'
subplot(1,3,2);
surf(log10(optimal_MLE_fit)); shading flat;
title 'MLE Poiss power law fit'
subplot(1,3,3);
residuals = maskeddata - optimal_MLE_fit;
residuals(maskeddata == -1) = nan;
surf(residuals); shading flat;
title 'MLE Poiss power law residuals'
colormap(flipud(pink));

%% Show the power law fit results.
figure;
image(org_circle_mask_data,'CDataMapping','Scaled');
colormap(pink);
hold on;
scatter(c_optimal(1),c_optimal(2),'y','filled');
scatter(c_logopt_optimal(1),c_logopt_optimal(2),'r','filled');
scatter(c_MLE_optimal(1),c_MLE_optimal(2),'g','filled');
legend('Natural space optimization','Log space optimization','Poisson MLE estimator');
title 'Uint16 image'

figure;
image(uint8(org_circle_mask_data),'CDataMapping','Scaled');
colormap(pink);
hold on;
scatter(c_optimal(1),c_optimal(2),'y','filled');
scatter(c_logopt_optimal(1),c_logopt_optimal(2),'r','filled');
scatter(c_MLE_optimal(1),c_MLE_optimal(2),'g','filled');
legend('Natural space optimization','Log space optimization','Poisson MLE estimator');
title 'Uint8 image'

%% (4) Bootstrap the Poisson MLE estimator 
% First check that the new function works
[linear_filtered_data] = makeFilteredLinearDataFromImageMatrix(maskeddata);
[origin_coords] = PoissonPowerlawBootfun(linear_filtered_data);%,c_optimal);
% scatter(origin_coords(1),origin_coords(2),'w','filled');
% Does seem to be working; I replicated the origin coordinates determined
% by the full Poisson method.
nboot = 30;
% Only run this when interested; takes up time.
%[ci,bootstat] = bootstrp(nboot,@PoissonPowerlawBootfun,linear_filtered_data);
% Compute [x,y,w,h] here (position vector);
% rectangle('Position',[x,y,w,h],'EdgeColor','g');

%% (5) Obtain beam position from COM of (a) graphene, and (b) hBN peaks 
% Use hexagonal window masks.
disp('Please specify hexagonal window masks for inner hBN.');
hexagonal_window_masks_hBN = getHexagonalWindow(data);
disp('Please specify hexagonal window masks for outer graphene.');
hexagonal_window_masks_graphene = getHexagonalWindow(data);

% Apply the convolution recipe 1.
kernel_radius = 5;
sigma = 0;
min_disk_distance = 8;
num_to_find_per_window = 1;
make_figures = 1;
[disk_locations, disk_intensities, com_coordinates_hBN] = braggDiskRegistrationRecipie_1(...
    data,hexagonal_window_masks_hBN,kernel_radius,sigma,min_disk_distance,num_to_find_per_window,make_figures);
num_to_find_per_window = 2;
[disk_locations, disk_intensities, com_coordinates_graphene] = braggDiskRegistrationRecipie_1(...
    data,hexagonal_window_masks_graphene,kernel_radius,sigma,min_disk_distance,num_to_find_per_window,make_figures);
% Note that subpixel resolution is not yet in effect here.



%% Compare global coordinates of all five methods.

[global_naturalopt_coords] = convertMaskedCoordsToOriginalCoords(c_optimal,Irng,Jrng)
[global_logopt_coords] = convertMaskedCoordsToOriginalCoords(c_logopt_optimal,Irng,Jrng)
[global_PoissonMLE_coords] = convertMaskedCoordsToOriginalCoords(c_MLE_optimal,Irng,Jrng)
global_hBN_coords = fliplr(com_coordinates_hBN)  % flipping needed because these come out as [ycoord,xcoord]
global_graphene_coords = fliplr(com_coordinates_graphene)

%% Global graphs
figure;
image(data,'CDataMapping','Scaled');
colormap(pink);
hold on;
scatter(global_naturalopt_coords(1),global_naturalopt_coords(2),'y','filled');
scatter(global_logopt_coords(1),global_logopt_coords(2),'r','filled');
scatter(global_PoissonMLE_coords(1),global_PoissonMLE_coords(2),'g','filled');
scatter(global_hBN_coords(1),global_hBN_coords(2),'m','filled');
scatter(global_graphene_coords(1),global_graphene_coords(2),'c','filled');
legend('Natural space optimization','Log space optimization','Poisson MLE estimator','Inner ring hBN center of mass','Outer ring graphene center of mass');
title 'Uint16 image'

figure;
image(uint8(data),'CDataMapping','Scaled');
colormap(pink);
hold on;
scatter(global_naturalopt_coords(1),global_naturalopt_coords(2),'y','filled');
scatter(global_logopt_coords(1),global_logopt_coords(2),'r','filled');
scatter(global_PoissonMLE_coords(1),global_PoissonMLE_coords(2),'g','filled');
scatter(global_hBN_coords(1),global_hBN_coords(2),'m','filled');
scatter(global_graphene_coords(1),global_graphene_coords(2),'c','filled');
legend('Natural space optimization','Log space optimization','Poisson MLE estimator','Inner ring hBN center of mass','Outer ring graphene center of mass');
title 'Uint8 image'

