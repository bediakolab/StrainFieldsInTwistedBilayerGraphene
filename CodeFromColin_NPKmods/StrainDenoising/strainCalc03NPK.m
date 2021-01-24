function [imageRGB,Emax,Emin,theta_piphase,Ediff,Edil] = strainCalc03NPK(exx,eyy,gxy)

% Colin Ophus - 2020 June
% Strain calculation from Nathanael's unwrapped displacement maps
%
% gxy here needs to be one half of the engineering strain!

MAKE_PLOTS = true;
tol = 1e-6;

% compute the principle strain directions, make color plots
% rangeE1 = [-0.05 0.1];
% rangeE2 = [0 -0.1];
rangeVal = [0.02 0.1];
rangeSat = [-0.02 0.02];
thetaShift = -0.2;

% Colin's original equations
% theta = atan(2*sStrain.Euv./(sStrain.Euu - sStrain.Evv)) / 2;
% E1 = (sStrain.Euu + sStrain.Evv)/2 ...
%     + sqrt((sStrain.Euu - sStrain.Evv).^2/4 + sStrain.Euv.^2);
% E2 = (sStrain.Euu + sStrain.Evv)/2 ...
%     - sqrt((sStrain.Euu - sStrain.Evv).^2/4 + sStrain.Euv.^2);
% NPK 06/14/2020 replaced these with the corrected values of the
% derivatives, now fed in via explicit arguments.
% Checked against the formulas here, and I believe my substitutions are
% correct: https://www.continuummechanics.org/principalstressesandstrains.html
theta = atan(2*gxy./(exx - eyy)) / 2;
Emax = (exx + eyy)/2 ...
    + sqrt((exx - eyy).^2/4 + gxy.^2);
Emin = (exx + eyy)/2 ...
    - sqrt((exx - eyy).^2/4 + gxy.^2);

% To obtain E1, E2, and theta with a range of pi rather than pi/2: cycle
% through with linear indices and figure out which rotation orientation has
% the larger component in the exx position.
Q = @(angle) [cos(angle), -sin(angle); sin(angle), cos(angle)];
theta_piphase = zeros(size(theta));
E1 = zeros(size(theta));
E2 = zeros(size(theta));
for i = 1:numel(theta)
    skipassert = false;
    thisexx = exx(i);
    thiseyy = eyy(i);
    thisgxy = gxy(i);
    thisangle = theta(i);
    thistensor = [thisexx,thisgxy;thisgxy,thiseyy];
    principaltensor = Q(thisangle)'*thistensor*Q(thisangle);
    principalrottensor90 = Q(thisangle+pi/2)'*thistensor*Q(thisangle+pi/2);
    if principaltensor(1,1) >= principaltensor(2,2)
        E1(i) = principaltensor(1,1);
        E2(i) = principaltensor(2,2);
        theta_piphase(i) = thisangle;
    elseif principalrottensor90(1,1) >= principalrottensor90(2,2)
        E1(i) = principalrottensor90(1,1);
        E2(i) = principalrottensor90(2,2);
        theta_piphase(i) = thisangle+pi/2;
    elseif isnan(principaltensor(1,1)) && isnan(principaltensor(2,2))
        E1(i) = nan;
        E2(i) = nan;
        theta_piphase(i) = nan;
        skipassert = true;
    else
        error('Tensor rotation has gone awry.');
    end
    % These give the same results, so the tensor rotation is working.
    % principalrottensor180 = Q(thisangle+pi)'*thistensor*Q(thisangle+pi);
    % principalrottensor270 = Q(thisangle+3*pi/2)'*thistensor*Q(thisangle+3*pi/2);
    % Assert that the rotation is in fact removing the shear strain
    % components
    if ~skipassert
        assert(abs(principaltensor(1,2))<tol);
        assert(abs(principaltensor(2,1))<tol);
        assert(abs(principalrottensor90(1,2))<tol);
        assert(abs(principalrottensor90(2,1))<tol);
    end
end

% If I'm understanding this correctly, this should be the same as the
% formulas Colin has above. The value in the above is simply to extend the
% phase of the angle.
% These always work, but screw up with the tear datasets.
% assert(nnz(~(abs(E1 - Emax) < tol)) == 0);
% assert(nnz(~(abs(E2 - Emin) < tol)) == 0);

% [min(theta(:)) max(theta(:))]*180/pi

% % Other combinations
% Colin's original equations
% Emax = max(abs(E1),abs(E2));
% Ediff = (E1 - E2);
% Edil = E1 + E2;
% % Rewritten with the new definitions
% Emax = max(abs(Emax),abs(Emin));
Ediff = (Emax - Emin);
Edil = Emax + Emin;

Ih = mod(theta/(pi/2) + thetaShift,1);

Iv = (-Emin - rangeVal(1)) / (rangeVal(2) - rangeVal(1));
Iv(:) = min(max(Iv,0),1);

Is = (Emax - rangeSat(1)) / (rangeSat(2) - rangeSat(1));
Is(:) = min(max(Is,0),1);

% Iv1 = (E1 - rangeE1(1)) / (rangeE1(2) - rangeE1(1));
% Iv2 = (E2 - rangeE2(1)) / (rangeE2(2) - rangeE2(1));
% Iv1(:) = min(max(Iv1,0),1);
% Iv2(:) = min(max(Iv2,0),1);

N = size(Emax);
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


if MAKE_PLOTS
    figure(1)
    clf
    imagesc([Emax Emin])
    axis equal off
    colormap(jet(256))
    colorbar
    caxis([-1 1]*0.05)
end    

figure
imagesc(imageRGB)
axis equal off
colormap(jet(256))
caxis([-1 1]*0.1)


end