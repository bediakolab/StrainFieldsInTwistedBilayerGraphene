function [sFit] = Gstrain03(sFit)

% Make vacuum probe template

threshMask = 0.02;
threshEdge = 0.005;

CBEDtemplate = sFit.CBEDvacuum;

% Shift position
sub = CBEDtemplate > (max(CBEDtemplate(:)) * threshMask);
[ya,xa] = meshgrid(1:sFit.stackSize(2),1:sFit.stackSize(1));
x0 = sum(sFit.CBEDvacuum(sub) .* xa(sub)) / sum(sFit.CBEDvacuum(sub));
y0 = sum(sFit.CBEDvacuum(sub) .* ya(sub)) / sum(sFit.CBEDvacuum(sub));
[qxa,qya] = makeFourierCoords(sFit.stackSize(1:2),1);
dxy = 1 - [x0 y0];
CBEDtemplate = ifft2(fft2(CBEDtemplate) ...
    .* exp(-2i*pi*(qxa*dxy(1) + qya*dxy(2))),'symmetric');
CBEDtemplate = max(CBEDtemplate,0);


% filter the edges of probe
sigma = 2;
k = fspecial('gaussian',2*ceil(3*sigma)+1,sigma);
CBEDfilt = convolve2(CBEDtemplate,k,'wrap');
s = max(CBEDfilt(:))*threshEdge;
maskEdge = 1 - exp(-CBEDfilt.^2/(2*s^2));
CBEDtemplate = CBEDtemplate .* maskEdge;


CBEDtemplateAmp = sqrt(max(CBEDtemplate,0));
scale = sqrt(max(CBEDtemplateAmp(:).^2));
CBEDtemplateAmp = CBEDtemplateAmp / scale;

sFit.CBEDtemplateAmp = CBEDtemplateAmp;


figure(11)
clf
imagesc(fftshift(CBEDtemplate))
hold on
hold off
axis equal off
colormap(jet(256))
set(gca,'position',[0 0 1 1])
% % caxis([0 10])

end