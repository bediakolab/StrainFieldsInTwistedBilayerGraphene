function [sStrain] = strainCalc02NPK(sStrain,strainDirRef)
% Colin Ophus - 2020 June
% Strain calculation from Nathanael's unwrapped displacement maps

% inputs
% strainDirRef - reference direction angle, in radians

% Remove this because we have set the strain range elsewhere.
% sStrain.strainRange = [-1 1]*0.1;
if strcmp(sStrain.strainRange,'auto')
    sStrain.rotationRange = 'auto';
else
    sStrain.rotationRange = [-1 1]*5 * pi/180;
end


if nargin < 2
    strainDirRef = 0;
end

% Masking
% maskStrain = ones(sStrain.sizeDisp(1:2)-2);
% NPK is not really sure what this is for
maskStrain = ones(size(sStrain.strainExx));

% inputs
Exx = sStrain.strainExx;
Eyy = sStrain.strainEyy;
% NPK changed the "Exy" to "Gxy" to avoid using the same symbol name twice.
Gxy = (sStrain.strainExy + sStrain.strainEyx) / 2;
sStrain.Erot = (sStrain.strainExy - sStrain.strainEyx)/2;

% Tensor rotation for strain
uP = [cos(strainDirRef) sin(strainDirRef)];
sStrain.Euu = uP(1)^2*Exx + 2*uP(1)*uP(2)*Gxy + uP(2)^2*Eyy;
sStrain.Evv = uP(1)^2*Eyy - 2*uP(1)*uP(2)*Gxy + uP(2)^2*Exx;
sStrain.Euv = (uP(1)^2-uP(2)^2)*Gxy - uP(1)*uP(2)*(Exx-Eyy);

% % color images
% EuuRGB = makeImage(sStrain.cMap,sStrain.Euu,sStrain.strainRange,maskStrain);
% EvvRGB = makeImage(sStrain.cMap,sStrain.Evv,sStrain.strainRange,maskStrain);
% EuvRGB = makeImage(sStrain.cMap,sStrain.Euv,sStrain.strainRange,maskStrain);
% phiRGB = makeImage(sStrain.cMap,sStrain.Erot,sStrain.rotationRange,maskStrain);
% NPK: make color images differently
figure
imagesc(sStrain.Euu*100); axis equal; set(gca,'yDir','normal'); cbh = colorbar; scaleColorMap(sStrain.cMap,0);
title('Exx rotated');
ylabel(cbh,'% strain');
hold on
line(24+[-1 1]*uP(2)*16,24+[-1 1]*uP(1)*16,'linewidth',2,'color','k');
figure
imagesc(sStrain.Evv); axis equal; set(gca,'yDir','normal'); cbh = colorbar; scaleColorMap(sStrain.cMap,0);
title('Eyy rotated');
ylabel(cbh,'% strain');
figure
imagesc(sStrain.Euv); axis equal; set(gca,'yDir','normal'); cbh = colorbar; scaleColorMap(sStrain.cMap,0);
title('Gxy rotated');
ylabel(cbh,'% strain');
figure
imagesc(rad2deg(sStrain.Erot)); axis equal; set(gca,'yDir','normal'); cbh = colorbar; scaleColorMap(sStrain.cMap,0);
title('Theta FBR rotated');
ylabel(cbh,'FBR (deg)');


% plotting
if false
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

end




function [imageRGB] = makeImage(cmap,sig,sigRange,mask)
sig = (sig - sigRange(1)) / (sigRange(2) - sigRange(1));
sig(:) = min(max(sig,0),1);
imageRGB = ind2rgb(round(255*sig)+1,cmap);
imageRGB(:) = rgb2hsv(imageRGB);
imageRGB(:,:,3) = mask;
imageRGB(:) = hsv2rgb(imageRGB);

end