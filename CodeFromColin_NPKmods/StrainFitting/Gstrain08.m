function [] = Gstrain08(sFit)



% Background
N = size(sFit.posRefine);
sig = sFit.posRefine(:,5,:,:);
I = reshape(permute(sig,[1 2 4 3]),[N(1)*N(3) N(4)])';
% Mean BG
I = squeeze(mean(sFit.posRefine(:,5,:,:),1));

% disks
sig = sFit.posRefine(:,6,:,:);
sig = sqrt(sFit.posRefine(:,7,:,:).^2 ...
    + sFit.posRefine(:,8,:,:).^2);
% sig = (atan2(sFit.posRefine(:,8,:,:),sFit.posRefine(:,7,:,:)));
% sig = atan(sFit.posRefine(:,8,:,:)./sFit.posRefine(:,7,:,:));
sig = abs((atan2(sFit.posRefine(:,8,:,:),sFit.posRefine(:,7,:,:))));

% I = reshape(permute(sig,[1 2 4 3]),[N(1)*N(3) N(4)])';
% I = squeeze(mean(sig,1));
% I = squeeze(median(sig,1));
% I = squeeze(sig(3,:,:,:));
% I = sqrt(squeeze( ...
%     median((sig-median(sig,1)).^2,1)));

mf = [1 1]*3;
I1 = medfilt2(squeeze(sig(1,:,:,:)));
I2 = medfilt2(squeeze(sig(2,:,:,:)));
I3 = medfilt2(squeeze(sig(3,:,:,:)));
I4 = medfilt2(squeeze(sig(4,:,:,:)));
I5 = medfilt2(squeeze(sig(5,:,:,:)));

% I = [ ...
%     zeros(40,10) ...
%     medfilt2(squeeze(sig(2,:,:,:))) ...
%     zeros(40,10) ...
%     medfilt2(squeeze(sig(3,:,:,:))) ...
%     zeros(40,10) ...
%     medfilt2(squeeze(sig(4,:,:,:))) ...
%     zeros(40,10) ...
%     medfilt2(squeeze(sig(5,:,:,:)))];

% NPK modification to get the dimensions worked out here.
I1 = I1';
I2 = I2';
I3 = I3';
I4 = I4';
I5 = I5';

I = [I1-median(I1(:)) ...
    zeros(40,5) ...
    I2-median(I2(:)) ...
    zeros(40,5) ...
    I3-median(I3(:)) ...
    zeros(40,5) ...
    I4-median(I4(:)) ...
    zeros(40,5) ...
    I5-median(I5(:))];

% I = medfilt2(I,[1 1]*3,'symmetric');

I = I - median(I(:));
I = I / sqrt(median(I(:).^2));

figure(12)
clf
imagesc(imresize(I,round(size(I).*[1 1.6]),'bilinear'))
axis equal off
% colormap(violetFire(256))
colormap(jet(256))
caxis([-1 1]*3)

end
