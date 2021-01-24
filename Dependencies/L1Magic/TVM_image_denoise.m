function [ denoised_image, rms_difference ] = TVM_image_denoise( noisy_image,...
    tolerance_fraction, plotflag, plotstring, colormaptype )
% Wrapper function to L1 Magic for total variation minimization denoising
% of an image (ROF framework, using identity matrix).
%
% Nathanael Kazmierczak, 04/03/2020

if (nargin < 3), plotflag = 0; end
if (nargin < 4), plotstring = 'image'; end
if (nargin < 5), colormaptype = gray; end

org_norm = norm(noisy_image(:));
xorg = noisy_image(:)/org_norm;
sorg_mean = mean(xorg);  % Not the true mean of the original image because we took the norm out first.
xorg = xorg - sorg_mean;
x0 = xorg;
b = xorg;
% Run in large mode with identity function handles
A = @(nv) nv;
At = @(kv) kv;

epsilon = norm(xorg)*tolerance_fraction;

[xp, ~] = tvqc_logbarrier(x0, A, At, b, epsilon);

% Sometimes the results have been going imaginary, so include this real
% call as a guard.
denoised_image_nc = reshape(real(xp),size(noisy_image));
denoised_image = denoised_image_nc + sorg_mean;
denoised_image = denoised_image * org_norm;
rms_difference = rms(rms(denoised_image - noisy_image));

if plotflag
    figure
    imagesc(noisy_image); set(gca, 'yDir','normal');
    colormap(colormaptype);
    colorbar
    title(sprintf('Original %s',plotstring)); axis equal;
    
    figure
    imagesc(denoised_image); set(gca, 'yDir','normal');
    colormap(colormaptype);
    colorbar
    title(sprintf('Denoised %s, error multiplier = %d',plotstring,tolerance_fraction));
    axis equal;
    
    fprintf('RMS image difference is %f.\n',rms_difference);
    fprintf('Epsilon quadratic constrain is %f.\n',epsilon);
end

end

