function [sStrain] = strainCalc01(unwrapped_disp_field)
tic
% Colin Ophus - 2020 June
% Strain calculation from Nathanael's unwrapped displacement maps
% Apply TVG denoising before calculating derivatives

% Inputs:
% unwrapped_disp_field - [Nx Ny 2] array containing displacements
%                        Assume displacement directions are x,y along dim 3
sStrain.flip_x_disp = true;
sStrain.flip_y_disp = false;
sStrain.flagNormalizeDisp = false;  % Remove mean plane from displacement maps

sStrain.strainRange = [-1 1]*0.10;  % for plotting

sStrain.TVG_padding = 20;
sStrain.TVG_iter = 64;  % number of iterations

sStrain.TVG_alpha = 1;  % smoothing parameter for TVG (smooth)
sStrain.TVG_beta = 0.5;  % smoothing parameter for TVG (TV)
sStrain.TVG_theta1 = 5;  % limiting parameter for TVG 
sStrain.TVG_theta2 = 5;  % limiting parameter for TVG 


sStrain.low_pass_sigma = 3;  % for comparison to the TVG denoising

% Strain plotting colormap
% cyan - white - red colormap
c = linspace(0,1,64)';
c0 = zeros(64,1);
c1 = ones(64,1);
sStrain.cMap = [ ...
    [c0 0.4+0.3*c c1];
    [c 0.7+0.3*c c1];
    [c1 1-c*0.4 1-c*0.4];
    [1-c*0.2 0.6-0.6*c 0.6-0.6*c];
    ];


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
sStrain.xDispTVG = TGV(sStrain.xDisp ,...
    sStrain.TVG_iter,....
    sStrain.TVG_padding,...
    sStrain.TVG_alpha,...
    sStrain.TVG_beta,...
    sStrain.TVG_theta1,...
    sStrain.TVG_theta2);
sStrain.yDispTVG = TGV(sStrain.yDisp ,...
    sStrain.TVG_iter,....
    sStrain.TVG_padding,...
    sStrain.TVG_alpha,...
    sStrain.TVG_beta,...
    sStrain.TVG_theta1,...
    sStrain.TVG_theta2);

% low pass filtering for comparison
k = fspecial('gaussian',...
    2*ceil(4*(sStrain.low_pass_sigma)),...
    sStrain.low_pass_sigma);
kNorm = 1./conv2(ones(sStrain.sizeDisp(1:2)),k,'same');
xDispLP = conv2(sStrain.xDisp,k,'same') .* kNorm;
yDispLP = conv2(sStrain.yDisp,k,'same') .* kNorm;


% derivatives - center difference
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

% plotting
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
imagesc([imagePlot ...
    zeros(2*sStrain.sizeDisp(1)-4,10)...
    imageCompare])
text(sStrain.sizeDisp(2)*0.5,10,'Exx',...
    'color','k','fontsize',16,'horizontalalign','center')
text(sStrain.sizeDisp(2)*1.5,10,'Exy',...
    'color','k','fontsize',16,'horizontalalign','center')
text(sStrain.sizeDisp(2)*0.5,10+sStrain.sizeDisp(1),'Eyx',...
    'color','k','fontsize',16,'horizontalalign','center')
text(sStrain.sizeDisp(2)*1.5,10+sStrain.sizeDisp(1),'Eyy',...
    'color','k','fontsize',16,'horizontalalign','center')

caxis(sStrain.strainRange)
axis equal off
colormap(sStrain.cMap)
set(gca,'position',[0 0 1 1])


toc
end