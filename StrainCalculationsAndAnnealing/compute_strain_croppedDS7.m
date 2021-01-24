% compute_strain_croppedDS7.m
%
% Script for exerting manual control to try to figure out why the
% displacements are not computing correctly.

addpath('/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/TEAM1Data/Dataset7RadialResidualsInterp');
load('03232020MagneticAnnealingExtendedZoneReconstruction.mat');
clear('figh1');
close all;


% annealing results.
figure
quiver(rgrid(:),cgrid(:),dfieldx(:),dfieldy(:),0);
hold on
quiver(fixed_emitters_indices(:,2),fixed_emitters_indices(:,1),fixed_emitters_directions(:,1),fixed_emitters_directions(:,2),0,'r');
set(gca,'YDir','normal');

figure
imagesc(dfieldx);
axis equal
title('Dfield x');
set(gca,'YDir','normal');
colormap jet
colorbar

figure
quiver(rgrid(:),cgrid(:),dfieldx(:),zeros(size(dfieldy(:))));
hold on
quiver(fixed_emitters_indices(:,2),fixed_emitters_indices(:,1),fixed_emitters_directions(:,1),fixed_emitters_directions(:,2),0,'r');
set(gca,'YDir','normal');
title('X component only');

figure
imagesc(dfieldy);
axis equal
title('Dfield y');
set(gca,'YDir','normal');
colormap jet
colorbar

%% Manually differentiate
exx = diff(dfieldx,1,2);
figure
imagesc(exx);
axis equal
title('exx');
set(gca,'YDir','normal');
colormap jet
colorbar
caxis([-0.5,0.5]);
title('raw exx');

exxfilt = diff(imgaussfilt(dfieldx,2),1,2);
figure
imagesc(exxfilt);
axis equal
title('exx');
set(gca,'YDir','normal');
colormap jet
colorbar
caxis([-0.5,0.5]);
title('gaussian exx');

%% What about just fitting the direction change?
angle = atan2(dfieldy,dfieldx);
[xchange,ychange] = gradient(angle);
changemag = sqrt(xchange.^2 + ychange.^2);
figure
imagesc(changemag);
axis equal
title('exx');
set(gca,'YDir','normal');
colormap jet
colorbar
caxis([0,0.25]);
title('angle change mag');

%% some sort of moving average approach
FILTER = 1;
AVERAGING_DOWNSAMPLE = 0;
int_factor = 2;

filtersize = 9;
angle = atan2(dfieldy,dfieldx);
% angle(angle > 1) = 0;
if FILTER
    movmean_angle = filter2(ones(filtersize)./49,angle);
    movmean_angle = trimArray(movmean_angle,4);
elseif AVERAGING_DOWNSAMPLE
    movmean_angle = averagingDownsample( angle, int_factor );
end
% movmean_angle = movmean(angle,5);
figure
imagesc(movmean_angle);
axis equal
title('exx');
set(gca,'YDir','normal');
colormap jet
colorbar
% caxis([0,0.25]);
title('Moving mean angle');

[xchange,ychange] = gradient(movmean_angle);
movemean_changemag = sqrt(xchange.^2 + ychange.^2);
figure
imagesc(movemean_changemag);
axis equal
title('exx');
set(gca,'YDir','normal');
colormap jet
colorbar
caxis([0,0.25]);
title('angle change mag for moving mean.');


%% Do the same moving average approach for the amplitude
amplitude = sqrt(dfieldy.^2 + dfieldx.^2);
% angle(angle > 1) = 0;
if FILTER
    movmean_amp = filter2(ones(filtersize)./49,amplitude);
    movmean_amp = trimArray(movmean_amp,4);
elseif AVERAGING_DOWNSAMPLE
    movmean_amp = averagingDownsample( amplitude, int_factor );
end
% movmean_angle = movmean(angle,5);
figure
imagesc(movmean_amp);
axis equal
title('exx');
set(gca,'YDir','normal');
colormap jet
colorbar
% caxis([0,0.25]);
title('Moving mean amplitude');

% This would need to be referenced against the local AB region -- converted
% back into reduced zone? But problem -- can we define a rotation and
% amplitude of the difference? Angle differences definitely don't work in
% extended zone, amplitude may not either (would have to work out the
% little triangles). Can we define amplitude and rotation in the
% equivalence class difference formalism? Seems like it should work -- this
% is how we know which AA region to reference against.


[xchange,ychange] = gradient(movmean_amp);
movemean_changemag = sqrt(xchange.^2 + ychange.^2);

figure
imagesc(movemean_changemag);
axis equal
title('exx');
set(gca,'YDir','normal');
colormap jet
colorbar
caxis([0,0.25]);
title('amplitude change mag for moving mean.');
