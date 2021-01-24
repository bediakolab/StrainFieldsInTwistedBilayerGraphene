% fit_AA_to_ellipse.m

addpath('/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/TEAM1Data/withMultistart/Dataset26');
load('DS26multistart_objectdata.mat');
addpath('/Volumes/NPKRes2');


displacement_field = m4.displacementContinuityFilter(0.3,9);
unfiltered_displacement_field = m4.DSC_fit_storage;
amplitude = sqrt(displacement_field(:,:,1).^2 + displacement_field(:,:,2).^2);
unfiltered_amplitude = sqrt(unfiltered_displacement_field(:,:,1).^2 + unfiltered_displacement_field(:,:,2).^2);
[ RGB_color_stack, HSV_color_stack ] = getCustomDisplacementColor( displacement_field, [], [], 1, 1 );
S = HSV_color_stack(:,:,2);
V = HSV_color_stack(:,:,3);
f = figure;
imagesc(amplitude);
title('Filtered amplitude');
colormap(gray);
colorbar;
set(gca,'yDir','normal');
% figure
% imagesc(unfiltered_amplitude);
% title('Unfiltered amplitude');
% colormap(gray);
% colorbar;
[figh] = m4.makeCustomDisplacementColorPlot([],[],0.3,9,1);
m4.makeCustomDisplacementColorPlot([],[],0.3,9,1);

thresh = 0.5;
amask = amplitude < thresh;
figure(figh);
hold on
h = imagesc(cat(3,~amask,~amask,zeros(size(amask))));
h.AlphaData = 0.5;

gauss_filter_thresh = 0.1;
edge_threshold = 10;
fAA = imgaussfilt(double(amask),5);
fAA(fAA < 0.1) = 0;
maxvals = imregionalmax(fAA);
lininds = find(maxvals);
[yv,xv] = ind2sub(size(fAA),lininds);
% ascertain proximity to array boundary
[r,c] = size(fAA);
adist = min([abs(1-yv),abs(r-yv),abs(1-xv),abs(c-xv)],[],2);
xv(adist < edge_threshold) = [];
yv(adist < edge_threshold) = [];
AA_centroids = horzcat(xv,yv);
nAA = size(AA_centroids,1);
xbase = 1:c;
ybase = 1:r;
[xspace,yspace] = meshgrid(xbase,ybase);
options.Display = 'iter';
scatter(AA_centroids(:,1),AA_centroids(:,2),100,'r','filled');

% Method 1: fit a circle to the thresholded region.
circlefit = false;
if circlefit
    predfun = @(x0,y0,r) isInCircle(xspace,yspace,x0,y0,r);
    fitstorage_circle = zeros(nAA,3);
    for i = 1:nAA
        thisAAcentroid = AA_centroids(i,:);
        thismask = amask;
        thismask(~isInCircle(xspace,yspace,thisAAcentroid(1),thisAAcentroid(2),50)) = 0;
        fitfun = @(params) nnz(thismask - predfun(params(1),params(2),params(3)));
        init_guesses = [thisAAcentroid,5];
        fitvals = fminsearch(fitfun,init_guesses,options);
        fitted_circle = predfun(fitvals(1),fitvals(2),fitvals(3));
        if false
            figure;
            imagesc(cat(3,thismask,zeros(size(thismask)),zeros(size(thismask))));
            hold on
            h = imagesc(cat(3,zeros(size(fitted_circle)),fitted_circle,zeros(size(fitted_circle))));
            h.AlphaData = 0.5;
            axis equal
        end
        fitstorage_circle(i,:) = fitvals;
    end % This does not appear to be terribly accurate -- it's OK
    
    figure(f)
    hold on
    for i = 1:nAA
        cent = fitstorage_circle(i,1:2);
%         cent2 = fitstorage_circle(i,2);
        r = fitstorage_circle(i,3);
        h=viscircles(cent,r,'Color','r');
    end
    axis equal
    
    % Get averaged AA diffraction patterns for Maddie
    for i = 1:nAA
        c = fitstorage_circle(i,:);
        this_mask = predfun(c(1),c(2),c(3));
        to_use = find(this_mask);
        to_use_subind = ind2sub(size(this_mask),to_use);
%         for j = 1:numel(to_use)
    end
    
    
    
end




% Method 3: fit a Gaussian to the actual amplitude data
% Still generate Gaussian over two dimensions
% Will fit in the linear data
gaussianfit = true;
if gaussianfit
    fitstorage_2DGauss = zeros(nAA,8);
    %     space = horzcat(xspace(:),yspace(:));
    %     gausspredfun = @(x0,y0,sigma1,sigma2,A,B,C) ellipticGaussianPredfun( space,x0,y0,sigma1,sigma2, A, B, C );
    
    % Set up the SP filter
    SP_hsv = [0,1,1;
        0.33,1,1;
        0.66,1,1];
    AB_hsv = [0,0,1;
        0.33,0,1;
        0.66,0,1];
    [ ~, HSV_color_stack ] = getCustomDisplacementColor( displacement_field, SP_hsv, AB_hsv, 1, 1 );
    H = HSV_color_stack(:,:,1);
    S = HSV_color_stack(:,:,2);
    V = HSV_color_stack(:,:,3);
        mask1 = S > 0.5;
        mask2 = V > 0.9;
        mask3 = mask1 & mask2;
    
    options = optimoptions('lsqnonlin');
    options.MaxFunEvals = 20000;
    options.MaxIter = 2000;
    options.CheckGradients = false;
    options.SpecifyObjectiveGradient = false;
    mask_radius = 25;
    clear c
    for i = 1:nAA
        i
        thisAAcentroid = AA_centroids(i,:);
        thismask = false(size(xspace));
        plotmask = false(size(xspace));
        thismask(isInCircle(xspace,yspace,thisAAcentroid(1),thisAAcentroid(2),mask_radius)) = true;
        thismask(mask3) = false;
        xrectmaskcoords = thisAAcentroid(1)-mask_radius:thisAAcentroid(1)+mask_radius;
        yrectmaskcoords = thisAAcentroid(2)-mask_radius:thisAAcentroid(2)+mask_radius;
        xrectmaskcoords(xrectmaskcoords<1) = [];
        yrectmaskcoords(yrectmaskcoords<1) = [];
        xrectmaskcoords(xrectmaskcoords>size(xspace,1)) = [];
        yrectmaskcoords(yrectmaskcoords>size(xspace,2)) = [];
        plotmask(yrectmaskcoords,xrectmaskcoords) = true;
%         residfun = @(c) ellipticGaussianResidfun( c, space, amplitude, thismask );%thismask(:).*(amplitude(:) - gausspredfun(c(1),c(2),c(3),c(4),c(5),c(6),c(7)));
        amplitude_to_fit = amplitude(thismask);
        xspace_touse = xspace(thismask);
        yspace_touse = yspace(thismask);
        space_tofit = horzcat(xspace_touse(:),yspace_touse(:));
        gausspredfun = @(c) ellipticGaussianPredfun2( space_tofit,c );
        residfun = @(c) ellipticGaussianResidfun2( c,space_tofit,amplitude_to_fit );
        rmsrfun = @(c) rms(residfun(c));
        
%         thisAAcentroid = fliplr(thisAAcentroid);
        init_guesses = [thisAAcentroid,4,4,0,1.4,1];
        lb = [thisAAcentroid-mask_radius,0,0,-1,0,0];
        ub = [thisAAcentroid+mask_radius,100,100,1,5,5];
        fitvals = lsqnonlin(residfun,init_guesses,lb,ub,options);
%         fitvals = fminsearch(rmsrfun,init_guesses,options);
        fitrmsr = rmsrfun(fitvals);
        
        
        amplitude_to_plot = amplitude(plotmask);
        amplitude_to_plug_in = amplitude(xrectmaskcoords,yrectmaskcoords);
        xspace_toplot = xspace(plotmask);
        yspace_toplot = yspace(plotmask);
        space_toplot = horzcat(xspace_toplot(:),yspace_toplot(:));
        fitted_gaussian_toplot = ellipticGaussianPredfun2( space_toplot,fitvals );
%         fitted_gaussian = gausspredfun(fitvals);
        fitted_gaussian_toplot = reshape(fitted_gaussian_toplot,[numel(yrectmaskcoords),numel(xrectmaskcoords)]);
        amplitude_to_plot = reshape(amplitude_to_plot,[numel(yrectmaskcoords),numel(xrectmaskcoords)]);
        if true
%             figure;
%             imagesc(cat(3,amplitude,zeros(size(amplitude_to_plot)),zeros(size(amplitude_to_plot))));
%             hold on
%             h = imagesc(cat(3,zeros(size(fitted_gaussian)),fitted_gaussian,zeros(size(fitted_gaussian))));
%             h.AlphaData = 0.5;
%             axis equal
            thismask_cropped = thismask(yrectmaskcoords,xrectmaskcoords);
            amplitude_to_plot2 = amplitude_to_plot;
            amplitude_to_plot2(~thismask_cropped) = nan;
            figure
            pcolor(xrectmaskcoords,yrectmaskcoords,-amplitude_to_plot2);
            axis equal
            shading flat
            title 'Amplitude data to fit'
            xlabel('xcoords');
            ylabel('ycoords');
            colorbar
            
            fitted_gaussian_toplot(~thismask_cropped) = nan;
            figure
            pcolor(xrectmaskcoords,yrectmaskcoords,-fitted_gaussian_toplot);
            axis equal
            shading flat
            title 'Gaussian fit'
            xlabel('xcoords');
            ylabel('ycoords');
            colorbar
            
            resid_toplot = ellipticGaussianResidfun2( fitvals,space_toplot,amplitude_to_plot(:) );
            resid_toplot = reshape(resid_toplot,size(thismask_cropped));
            resid_toplot(~thismask_cropped) = nan;
            figure
            pcolor(xrectmaskcoords,yrectmaskcoords,resid_toplot);
            axis equal
            shading flat
            title 'Residuals'
            xlabel('xcoords');
            ylabel('ycoords');
            colorbar
        end
        fitstorage_2DGauss(i,:) = [fitvals,fitrmsr];
    end
    fitstorage_2DGauss
    
    % Plot
%     figh1 = m4.makeCustomDisplacementColorPlot([],[],0.3,9,1);


    


figure(f)
    hold on
    angle_axis_storage = zeros(nAA,3);
    for i = 1:nAA
        % fsolve to back out semimajor/minor axes and rotation angle
        FWHM1 = fitstorage_2DGauss(i,3);
        FWHM2 = fitstorage_2DGauss(i,4);
        p = fitstorage_2DGauss(i,5);
        s1 = FWHM1/(2*sqrt(2*log(2)));
        s2 = FWHM2/(2*sqrt(2*log(2)));
        
        k = 2*log(2)*(1-p^2);
        eqn1 = @(alpha,a,b) cos(alpha)^2/(a^2) + sin(alpha)^2/(b^2) - 1/(k*s1^2);
        eqn2 = @(alpha,a,b) sin(2*alpha)*(1/(a^2) - 1/(b^2)) + 2*p/(k*s1*s2);
        eqn3 = @(alpha,a,b) sin(alpha)^2/(a^2) + cos(alpha)^2/(b^2) - 1/(k*s2^2);
        eqns = @(x) [eqn1(x(1),x(2),x(3)),eqn2(x(1),x(2),x(3)),eqn3(x(1),x(2),x(3))];
        x0 = [p,FWHM1,FWHM2];
        options = optimset('lsqnonlin');
        options.TolFun = 1e-10;
        % set bounds because there seems to be some sort of ambiguity left
        % in these equations
        if p > 0
            lb = [0,0,0];
            ub = [pi/2,inf,inf];
        else
            lb = [-pi/2,0,0];
            ub = [0,inf,inf];
        end
        solvevals = lsqnonlin(eqns,x0,lb,ub,options)
        residuals = eqns(solvevals)
        assert(rms(residuals) < 1e-5);
        
        semimajor = max(solvevals(2:3));
        semiminor = min(solvevals(2:3));
        angle = solvevals(1);  % no attempt here to phase wrap, as we don't really care about the angle value anyways.
        angle_axis_storage(i,:) = [angle,semimajor,semiminor];
%         ang = fitstorage_2DGauss(i,7);
%         ra = fitstorage_2DGauss(i,3);
%         rb = fitstorage_2DGauss(i,4);
        % these are correct
        x0 = fitstorage_2DGauss(i,1);
        y0 = fitstorage_2DGauss(i,2);
        h=ellipse(semimajor,semiminor,angle,x0,y0);
    end
    axis equal
end

