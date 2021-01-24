function [imageRGB,E1,E2,theta] = strainCalc03(sStrain)

% Colin Ophus - 2020 June
% Strain calculation from Nathanael's unwrapped displacement maps

% compute the principle strain directions, make color plots
% rangeE1 = [-0.05 0.1];
% rangeE2 = [0 -0.1];
rangeVal = [0.02 0.1];
rangeSat = [-0.02 0.02];
thetaShift = -0.2;


theta = atan(2*sStrain.Euv./(sStrain.Euu - sStrain.Evv)) / 2;
E1 = (sStrain.Euu + sStrain.Evv)/2 ...
    + sqrt((sStrain.Euu - sStrain.Evv).^2/4 + sStrain.Euv.^2);
E2 = (sStrain.Euu + sStrain.Evv)/2 ...
    - sqrt((sStrain.Euu - sStrain.Evv).^2/4 + sStrain.Euv.^2);

% [min(theta(:)) max(theta(:))]*180/pi

% % Other combinations
Emax = max(abs(E1),abs(E2));
Ediff = (E1 - E2);
Edil = E1 + E2;

Ih = mod(theta/(pi/2) + thetaShift,1);

Iv = (-E2 - rangeVal(1)) / (rangeVal(2) - rangeVal(1));
Iv(:) = min(max(Iv,0),1);

Is = (E1 - rangeSat(1)) / (rangeSat(2) - rangeSat(1));
Is(:) = min(max(Is,0),1);

% Iv1 = (E1 - rangeE1(1)) / (rangeE1(2) - rangeE1(1));
% Iv2 = (E2 - rangeE2(1)) / (rangeE2(2) - rangeE2(1));
% Iv1(:) = min(max(Iv1,0),1);
% Iv2(:) = min(max(Iv2,0),1);

N = size(E1);
imageRGB = ones(N(1),N(2),3);
imageRGB(:,:,1) = Ih;
imageRGB(:,:,2) = Is;
imageRGB(:,:,3) = Iv;
imageRGB(:) = hsv2rgb(imageRGB);

% imageRGB_E1 = ones(N(1),N(2),3);
% imageRGB_E2 = ones(N(1),N(2),3);
% imageRGB_E1(:,:,1) = Ih;
% imageRGB_E2(:,:,1) = Ih;
% imageRGB_E1(:,:,3) = Iv1;
% imageRGB_E2(:,:,3) = Iv2;
% 
% imageRGB_E1(:) = hsv2rgb(imageRGB_E1);
% imageRGB_E2(:) = hsv2rgb(imageRGB_E2);



figure(1)
clf
imagesc([E1 E2])
axis equal off
colormap(jet(256))
colorbar
caxis([-1 1]*0.05)


figure(2)
clf
imagesc(imageRGB)
axis equal off
colormap(jet(256))
caxis([-1 1]*0.1)


end