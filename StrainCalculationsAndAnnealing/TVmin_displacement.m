% TVmin_displacement.m
%
% 04/03/2020

% load in an extended zone annealing reconstruction in dfield

dfieldx = dfield(:,:,1);
dfieldy = dfield(:,:,2);
xdisp = dfieldx(:)/norm(dfieldx(:));
xorg = xdisp - mean(xdisp);
x0 = xorg;
b = xorg;
% x = reshape(I,N,1);
A = @(nv) nv;
At = @(kv) kv;

frac = 0.2;
epsilon = norm(xorg)*frac;

[xp, tp] = tvqc_logbarrier(x0, A, At, b, epsilon);

xp = reshape(real(xp),size(dfieldx));
figure
imagesc(xp); set(gca, 'yDir','normal');
colorbar
title(sprintf('Denoised x displacement field, error multiplier = %d',frac));
figure
xorg2 = reshape(b,size(dfieldx));
imagesc(xorg2); set(gca, 'yDir','normal');
colorbar
title('Original x displacement field');

rms(rms(xp - xorg2))
epsilon



[exx,exy] = gradient(xp);
figure;
imagesc(exx); set(gca, 'yDir','normal'); colorbar;
title('exx');
figure
imagesc(exy); set(gca, 'yDir','normal'); colorbar;
title('exy');

