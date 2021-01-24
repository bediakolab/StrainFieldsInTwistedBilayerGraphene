function [sStrain] = strainCalc01NPK(unwrapped_disp_field,sStrain,make_plots)
tic
% Colin Ophus - 2020 June
% Strain calculation from Nathanael's unwrapped displacement maps
% Apply TVG denoising before calculating derivatives
%
% Modified by Nathanael Kazmerczak 2020 June.
% sStrain options will now be set in the caller function, which is
% FourDSTEM_Analysis_Engine.computeStrainMaps3().

% Inputs:
% unwrapped_disp_field - [Nx Ny 2] array containing displacements
%                        Assume displacement directions are x,y along dim 3




% Strain plotting colormap
% cyan - white - red colormap
% c = linspace(0,1,64)';
% c0 = zeros(64,1);
% c1 = ones(64,1);
% sStrain.cMap = [ ...
%     [c0 0.4+0.3*c c1];
%     [c 0.7+0.3*c c1];
%     [c1 1-c*0.4 1-c*0.4];
%     [1-c*0.2 0.6-0.6*c 0.6-0.6*c];
%     ];
% NPK: use the Jon King colormap
load('/Users/nathanaelkazmierczak/Dropbox/SharedOS/NPKHardDriveBackupFall2017/workspace/BediakoResearch/matlab4DSTEM/Dependencies/DivergingColormap/JonKing93-scaleColorMap-23ed82e/Demo/demo_colormap.mat');
sStrain.cMap = cmap;       

% input data
sStrain.xDisp = unwrapped_disp_field(:,:,1);
sStrain.yDisp = unwrapped_disp_field(:,:,2);
sStrain.sizeDisp = size(unwrapped_disp_field);

% Flipping
if sStrain.flip_x_disp == true
    sStrain.xDisp(:) = -sStrain.xDisp;
end
if sStrain.flip_y_disp == true
    sStrain.yDisp(:) = -sStrain.yDisp;
end

% normalization
if sStrain.flagNormalizeDisp == true
    % subtract plane
    
    [ya,xa] = meshgrid(1:sStrain.sizeDisp(2),1:sStrain.sizeDisp(1));
    basis = [ones(prod(sStrain.sizeDisp(1:2)),1) xa(:) ya(:)];
    
    xBeta = basis \ sStrain.xDisp(:);
    yBeta = basis \ sStrain.yDisp(:);
    
    sStrain.xDisp = sStrain.xDisp ....
        - reshape(basis*xBeta,sStrain.sizeDisp(1:2));
    sStrain.yDisp = sStrain.yDisp ....
        - reshape(basis*yBeta,sStrain.sizeDisp(1:2));
end

% run TGV denoising
if sStrain.TVG_iter > 0
disp('Beginning TGV denoising...');
sStrain.xDispTVG = TGV_NPK(sStrain.xDisp ,...
    sStrain.TVG_iter,....
    sStrain.TVG_padding,...
    sStrain.TVG_alpha,...
    sStrain.TVG_beta,...
    sStrain.TVG_theta1,...
    sStrain.TVG_theta2);
sStrain.yDispTVG = TGV_NPK(sStrain.yDisp ,...
    sStrain.TVG_iter,....
    sStrain.TVG_padding,...
    sStrain.TVG_alpha,...
    sStrain.TVG_beta,...
    sStrain.TVG_theta1,...
    sStrain.TVG_theta2);
disp('TGV denoising complete.');
else
    sStrain.xDispTVG = sStrain.xDisp;
    sStrain.yDispTVG = sStrain.yDisp;
end

% low pass filtering for comparison
k = fspecial('gaussian',...
    2*ceil(4*(sStrain.low_pass_sigma)),...
    sStrain.low_pass_sigma);
kNorm = 1./conv2(ones(sStrain.sizeDisp(1:2)),k,'same');
xDispLP = conv2(sStrain.xDisp,k,'same') .* kNorm;
yDispLP = conv2(sStrain.yDisp,k,'same') .* kNorm;


% derivatives - center difference
% NPK note: Colin originally has the centered finite difference plugged in
% manually, so I'm not sure that it respects the stepsize.
%
% NOTE: the gradient function computes derivatives in the way I'm used to
% thinking, with the horizontal direction treated as "x"
USE_GRADIENT_FUNCTION = true;
if USE_GRADIENT_FUNCTION
    % The way NPK originally computes the derivatives.
    uX = sStrain.xDispTVG;
    uY = sStrain.yDispTVG;
    [dUinterX_dorgx,dUinterX_dorgy] = gradient(uX,sStrain.stepsize);
    [dUinterY_dorgx,dUinterY_dorgy] = gradient(uY,sStrain.stepsize);
    sStrain.strainExx = dUinterX_dorgx;
    sStrain.strainExy = dUinterX_dorgy;
    sStrain.strainEyx = dUinterY_dorgx;
    sStrain.strainEyy = dUinterY_dorgy;
    [ExxLP,ExyLP] = gradient(xDispLP,sStrain.stepsize);
    [EyxLP,EyyLP] = gradient(yDispLP,sStrain.stepsize);
else
    vx = 2:(sStrain.sizeDisp(1)-1);
    vy = 2:(sStrain.sizeDisp(2)-1);
    sStrain.strainExx = (sStrain.xDispTVG(vx+1,vy) - sStrain.xDispTVG(vx-1,vy))/2;
    sStrain.strainExy = (sStrain.xDispTVG(vx,vy+1) - sStrain.xDispTVG(vx,vy-1))/2;
    sStrain.strainEyx = (sStrain.yDispTVG(vx+1,vy) - sStrain.yDispTVG(vx-1,vy))/2;
    sStrain.strainEyy = (sStrain.yDispTVG(vx,vy+1) - sStrain.yDispTVG(vx,vy-1))/2;
    
    ExxLP = (xDispLP(vx+1,vy) - xDispLP(vx-1,vy))/2;
    ExyLP = (xDispLP(vx,vy+1) - xDispLP(vx,vy-1))/2;
    EyxLP = (yDispLP(vx+1,vy) - yDispLP(vx-1,vy))/2;
    EyyLP = (yDispLP(vx,vy+1) - yDispLP(vx,vy-1))/2;
    step_correct = true;
    if step_correct
        % Divide by the stepsize, because that is the correction to the
        % unit step that Colin has been assuming.
        sStrain.strainExx = sStrain.strainExx./sStrain.stepsize;
        sStrain.strainExy = sStrain.strainExy./sStrain.stepsize;
        sStrain.strainEyx = sStrain.strainEyx./sStrain.stepsize;
        sStrain.strainEyy = sStrain.strainEyy./sStrain.stepsize;
        ExxLP = ExxLP./sStrain.stepsize;
        ExyLP = ExyLP./sStrain.stepsize;
        EyxLP = EyxLP./sStrain.stepsize;
        EyyLP = EyyLP./sStrain.stepsize;
    end
end



% plotting
if make_plots
    imagePlot = [ ...
        sStrain.strainExx sStrain.strainExy;
        sStrain.strainEyx sStrain.strainEyy;
        ];
    imageCompare = [ ...
        ExxLP ExyLP;
        EyxLP EyyLP;
        ];
    
    
    
    figure(21)
    clf
    if USE_GRADIENT_FUNCTION
        % This doesn't lose edge pixels like Colin's manual derivatives do.
        % Therefore take out the "-4" in the size.
        imagesc([imagePlot ...
            zeros(2*sStrain.sizeDisp(1),10)...
            imageCompare]*100);  % NPK added the factor of 100 so that this will be a true strain % (it is nm/nm * 100)
    else
        imagesc([imagePlot ...
            zeros(2*sStrain.sizeDisp(1)-4,10)...
            imageCompare]*100);  % NPK added the factor of 100 so that this will be a true strain % (it is nm/nm * 100)
    end
    text(sStrain.sizeDisp(2)*0.5,10,'Exx',...
        'color','k','fontsize',16,'horizontalalign','center')
    text(sStrain.sizeDisp(2)*1.5,10,'Exy',...
        'color','k','fontsize',16,'horizontalalign','center')
    text(sStrain.sizeDisp(2)*0.5,10+sStrain.sizeDisp(1),'Eyx',...
        'color','k','fontsize',16,'horizontalalign','center')
    text(sStrain.sizeDisp(2)*1.5,10+sStrain.sizeDisp(1),'Eyy',...
        'color','k','fontsize',16,'horizontalalign','center')
    
    % NPK: add some legends here
    text(sStrain.sizeDisp(2)*1,-15,'Total Generalized Variation (TGV) strain maps',...
        'color','k','fontsize',16,'horizontalalign','center')
    text(sStrain.sizeDisp(2)*3,-15,'Gaussian lowpass filter strain maps',...
        'color','k','fontsize',16,'horizontalalign','center')
    
    if ~strcmp(sStrain.strainRange,'auto')
        caxis(sStrain.strainRange)
        colormap(sStrain.cMap)
    else
        load('/Users/nathanaelkazmierczak/Dropbox/SharedOS/NPKHardDriveBackupFall2017/workspace/BediakoResearch/matlab4DSTEM/Dependencies/DivergingColormap/JonKing93-scaleColorMap-23ed82e/Demo/demo_colormap.mat');
        scaleColorMap( cmap, 0 );
        sStrain.cMap = cmap;
    end
    axis equal off
    
    set(gca,'position',[0 0 1 1])
    % NPK add colorbar
    cbh = colorbar;
    ylabel(cbh,'Strain %');
    
end

toc
end