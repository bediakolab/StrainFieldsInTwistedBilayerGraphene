function [sStrain] = strainCalc02(sStrain,strainDirRef)
% Colin Ophus - 2020 June
% Strain calculation from Nathanael's unwrapped displacement maps

% inputs
% strainDirRef - reference direction angle, in radians

sStrain.strainRange = [-1 1]*0.1;
sStrain.rotationRange = [-1 1]*5 * pi/180;


if nargin < 2
    strainDirRef = 0;
end

% Masking
maskStrain = ones(sStrain.sizeDisp(1:2)-2);

% inputs
Exx = sStrain.strainExx;
Eyy = sStrain.strainEyy;
Exy = (sStrain.strainExy + sStrain.strainEyx) / 2;
sStrain.Erot = (sStrain.strainExy - sStrain.strainEyx)/2;

% Tensor rotation for strain
uP = [cos(strainDirRef) sin(strainDirRef)];
sStrain.Euu = uP(1)^2*Exx + 2*uP(1)*uP(2)*Exy + uP(2)^2*Eyy;
sStrain.Evv = uP(1)^2*Eyy - 2*uP(1)*uP(2)*Exy + uP(2)^2*Exx;
sStrain.Euv = (uP(1)^2-uP(2)^2)*Exy - uP(1)*uP(2)*(Exx-Eyy);

% color images
EuuRGB = makeImage(sStrain.cMap,sStrain.Euu,sStrain.strainRange,maskStrain);
EvvRGB = makeImage(sStrain.cMap,sStrain.Evv,sStrain.strainRange,maskStrain);
EuvRGB = makeImage(sStrain.cMap,sStrain.Euv,sStrain.strainRange,maskStrain);
phiRGB = makeImage(sStrain.cMap,sStrain.Erot,sStrain.rotationRange,maskStrain);


% plotting
imagePlot = [ ...
    EuuRGB EvvRGB;
    EuvRGB phiRGB];

figure(21)
clf
imagesc(imagePlot)
text(sStrain.sizeDisp(2)*0.5,10,'Euu',...
    'color','k','fontsize',16,'horizontalalign','center')
text(sStrain.sizeDisp(2)*1.5,10,'Evv',...
    'color','k','fontsize',16,'horizontalalign','center')
text(sStrain.sizeDisp(2)*0.5,10+sStrain.sizeDisp(1),'Euv',...
    'color','k','fontsize',16,'horizontalalign','center')
text(sStrain.sizeDisp(2)*1.5,10+sStrain.sizeDisp(1),'Phi',...
    'color','k','fontsize',16,'horizontalalign','center')

line(24+[-1 1]*uP(2)*16,24+[-1 1]*uP(1)*16,'linewidth',2,'color','k')


% caxis(sStrain.strainRange)
axis equal off
colormap(sStrain.cMap)
set(gca,'position',[0 0 1 1])

end




function [imageRGB] = makeImage(cmap,sig,sigRange,mask)
sig = (sig - sigRange(1)) / (sigRange(2) - sigRange(1));
sig(:) = min(max(sig,0),1);
imageRGB = ind2rgb(round(255*sig)+1,cmap);
imageRGB(:) = rgb2hsv(imageRGB);
imageRGB(:,:,3) = mask;
imageRGB(:) = hsv2rgb(imageRGB);

end