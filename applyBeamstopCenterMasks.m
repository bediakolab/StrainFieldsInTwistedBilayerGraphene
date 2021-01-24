function [maskeddata] = applyBeamstopCenterMasks(data,masks)
% Function that uses masks produced by makeBeamstopCenterMasks.m to mask an
% image for power law fitting of the center of the image.
%
% -1 is the masking signal, as interpreted by the Poisson MLE functions,
% etc.
%
% Nathanael Kazmierczak, 08/11/2019
% Bediako Lab, UC Berkeley

THRESHOLD = 50;

maskeddata = double(data);
maskeddata(~masks.circle_mask) = -1;
maskeddata = maskeddata(masks.Irng(1):masks.Irng(2),masks.Jrng(1):masks.Jrng(2));
maskeddata(masks.beamstop_mask & maskeddata < THRESHOLD) = -1;
maskeddata(masks.inner_circle_mask & masks.inner_beamstop_mask) = -1;
if isfield(masks,'teeny_circle_mask')
    maskeddata(teeny_circle_mask) = -1;
end

end

