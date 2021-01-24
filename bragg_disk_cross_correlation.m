% bragg_disk_cross_correlation.m

% Attempting to do in Matlab what has not been working in Python.
% First attempts
% Nathanael Kazmierczak

addpath /Volumes/Lexar
info = h5info('3__70x70_ss4nm_0p05s_spot9_alpha=0p55_bin=4_CL=130_80kV_C2=40um__cropped.h5');
h5disp('3__70x70_ss4nm_0p05s_spot9_alpha=0p55_bin=4_CL=130_80kV_C2=40um__cropped.h5','/4DSTEM_experiment/data/datacubes/datacube_0/data');

% Attempt to read in exactly three diffraction patterns
start = [1,1,1,1];
count = [512,512,3,1];
data = h5read('3__70x70_ss4nm_0p05s_spot9_alpha=0p55_bin=4_CL=130_80kV_C2=40um__cropped.h5','/4DSTEM_experiment/data/datacubes/datacube_0/data',start,count);
% imshow(uint8(data(:,:,1)));
% figure
% contourf(data(:,:,1),200);

figure;
pcolor(uint8(data(:,:,1))); shading flat;

DP = data(:,:,1); % Diffraction pattern of interest.

%% Generate synthetic probe kernel
idx = round(512/2);
r = 5;
xbase = 1:512;
ybase = 1:512;
[xspace,yspace] = meshgrid(xbase,ybase);
probe_kernel = isInCircle(xspace,yspace,idx,idx,r);
figure
pcolor(probe_kernel); shading flat;

%% Cross-correlate with the diffraction pattern
% note that here the images are the same size.

corr = xcorr2(double(uint8(DP)),double(probe_kernel));
figure; image(corr); shading flat; title 'Cross correlation matrix';

% Need to map back onto the original image. This line deserves further
% scrutiny in the future to ascertain whether it is correct.
overlay = corr(255:766,255:766);
figure; image(uint8(DP),'CDataMapping','scaled'); title 'Overlay'; 
%ax = gca;
hold all; contour(overlay/100);
colormap jet;


BW = imregionalmax(overlay);
figure;
pcolor(BW); shading flat; title 'Found peaks raw'
sigma = 1;
gauss_smooth = imgaussfilt(overlay,1);
line_profile = improfile(overlay,[256,0],[256,0]);
figure;
plot(line_profile);
coeff_guesses = [1e10,-4];

[fit1,optimalparams1] = getPowerLawFit(line_profile',1:numel(line_profile),coeff_guesses);
figure;
plot(1:numel(line_profile),line_profile,'ob',1:numel(line_profile),fit1,'-r');

%% Now use the background correction radially:
radii = sqrt((xspace - 256).^2 + (yspace - 256).^2);
powerlaw_pred = @(c,x) c(1).*(x).^(c(2));
baseline = powerlaw_pred(optimalparams1,radii);
baseline(baseline > 10000) = 10000; % because otherwise it goes super negative near (0,0).
figure;
surf(baseline); shading flat; title 'Power law baseline'
baselined_corr_data = overlay - baseline;
figure;
surf(baselined_corr_data); shading flat; colormap hsv; title 'Baselined correlation data'

% Decree also that any pixels less than 50 distance from (270,260) should be set to zero.
% This is closer to where I visually think the middle is.
radii2 = sqrt((xspace - 270).^2 + (yspace - 260).^2);
baselined_corr_data_annular = baselined_corr_data;
baselined_corr_data_annular(radii2 < 80) = 0;
figure;
surf(baselined_corr_data_annular); shading flat; colormap hsv; title 'Baselined, annularized correlation data'


gauss_smooth2 = imgaussfilt(baselined_corr_data,0.1);
BW2 = imregionalmax(gauss_smooth2);
BW2(baselined_corr_data_annular < 1000) = 0;
figure;
spy(flipud(BW2)); shading flat; title 'Found peaks processed'


% Implement annular filters
x0 = 270;
y0 = 260;
rinner = 80;
router = 100;
[annulus_mask_1] = isInAnnulus(xspace,yspace,x0,y0,rinner,router);
rinner2 = 145;
router2 = 170;
[annulus_mask_2] = isInAnnulus(xspace,yspace,x0,y0,rinner2,router2);
baselined_corr_data_annular2 = baselined_corr_data_annular;
baselined_corr_data_annular2(~annulus_mask_2 & ~annulus_mask_1) = 0;
figure;
contour(baselined_corr_data_annular2,100); shading flat; colormap parula; title 'Baselined, more annulus correlation data'

%% Can do hexagonal window masking
[hexagonal_window_masks] = getHexagonalWindow(overlay);
masked_data = baselined_corr_data;
masked_data(~sum(hexagonal_window_masks,3)) = 0;
figure; surf(masked_data); shading flat;

%% Not using the power law fit...
sigma = 0.1;
smootheddata = imgaussfilt(overlay,sigma);
BW3 = imregionalmax(smootheddata);
masked_maxima = BW3;
masked_maxima(~sum(hexagonal_window_masks,3)) = 0;
lininds_maxima = find(masked_maxima);
[xinds,yinds] = ind2sub(size(masked_maxima),lininds_maxima);
intensities = overlay(lininds_maxima);

combodata = [intensities,xinds,yinds];
sortedcombodata = sortrows(combodata);

figure; 
plotmask = sum(hexagonal_window_masks,3);
plotmask(plotmask == 0) = 0.1;
plotdata = plotmask .* log10(smootheddata);
pcolor(plotdata); shading flat;
ax = gca;
hold on;
h = scatter(ax,yinds,xinds,'filled'); colormap(ax,pink);  % apparently needed to get them to line up correctly

% Now, find the two most intense peaks in each window.
maximastorage = zeros(12,3);
for i = 1:6
    this_mask = hexagonal_window_masks(:,:,i);
    these_maxima = BW3 & this_mask;
    lininds_maxima = find(these_maxima);
    [xinds,yinds] = ind2sub(size(these_maxima),lininds_maxima);
    intensities = baselined_corr_data(lininds_maxima);
    combodata = [intensities,xinds,yinds];
    sortedcombodata = flipud(sortrows(combodata));
    
    maximastorage(i*2-1:i*2,:) = sortedcombodata(1:2,:);
end

% Plot on top of the maxima map from correlation and smoothing
figure; 
plotmask = sum(hexagonal_window_masks,3);
plotmask(plotmask == 0) = 0.1;
plotdata = plotmask .* log10(smootheddata);
pcolor(plotdata); shading flat;
ax = gca;
hold on;
h = scatter(ax,maximastorage(:,3),maximastorage(:,2),'filled'); colormap(ax,pink);

figure; 
plotmask = sum(hexagonal_window_masks,3);
plotmask(plotmask == 0) = 0.1;
plotdata = plotmask .* double(uint8(DP));
pcolor(plotdata); shading flat;
ax = gca;
hold on;
h = scatter(ax,maximastorage(:,3),maximastorage(:,2),'filled'); colormap(ax,pink);






% % % h2 = pcolor(plotdata); shading flat;
% % % hold on;
% % % h1 = pcolor(smootheddata); shading flat;
% % % set(h1,'facealpha',0);
% % % ax = gca;
% % % h1.AlphaData = sum(hexagonal_window_masks,3);
% % % set(h1,'facealpha',0.5);
% % % hold on;
% % % h2 = pcolor(ax,sum(hexagonal_window_masks,3)); shading flat;
% set(h2,'facealpha',0.5);


% figure; scatter(masked_maxima); shading flat;
% % % % 
% % % % figure; pcolor(masked_maxima); shading flat;
% % % % figure; scatter(sortedcombodata(1:10,2),sortedcombodata(1:10,3)); shading flat;
% % 
% % start = [1,1,1,1];
% % count = [512,512,20,1];
% % data = h5read('3__70x70_ss4nm_0p05s_spot9_alpha=0p55_bin=4_CL=130_80kV_C2=40um__cropped.h5','/4DSTEM_experiment/data/datacubes/datacube_0/data',start,count);


