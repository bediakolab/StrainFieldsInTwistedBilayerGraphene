% trigFittingFunction_testing_script.m
%
% It seems there is something subtley wrong with the bounds. This is to try
% to figure out why.
%
% Nathanael Kazmierczak, 02/02/2020

filename = '7_80x80_ss=2nm_C2=40um_alpha=1mrad_spot7_200ms_CL=130 bin=4_60kV_NPKnormalized_RSbin2.h5';
cd ..
addpath(genpath(pwd));
cd DSC_BlinkingFit_Code
addpath('/Volumes/NPKResearch/BediakoLab/Maddie4DSTEM/20200123');
data = h5read(filename,'/4DSTEM_experiment/data/datacubes/datacube_0/data',[1,1,20,20],[512,512,1,1]);
plotDP(data);
            
loadpath = '/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/DSC_BlinkingFit/Results/0123_dataset7_RSbin2_try2padremoved';
thisdir = pwd;
cd(loadpath);
load('objectData.mat');
cd(thisdir);

use_teeny_mask = 0;
teeny_mask_radius_factor = 0;
draw_hard_delete_region = 2;
THRESHOLD = 150;
[beamstop_masks,Irng,Jrng,maskeddata] = makeBeamstopCenterMasks(data,use_teeny_mask,teeny_mask_radius_factor,draw_hard_delete_region,THRESHOLD);

powerLawFit = @(c) getFullRadialPowerLawFit(c,maskeddata);
y0_init = round(mean([1,size(maskeddata,1)]));
x0_init = round(mean([1,size(maskeddata,2)]));
log10A_init = 15;
B_init = -4;
c_initial_guesses = [x0_init,y0_init,log10A_init,B_init];
options = optimset;
options.Display = 'iter';
options.TolFun = 1e-7;
options.TolX = 1e-7;

c_optimal = fminsearch(powerLawFit,c_initial_guesses,options);

optimal_fit = getFullRadialPowerLawPred(c_optimal,size(maskeddata));
figure; subplot(1,3,1);
nan_masked_data = maskeddata;
nan_masked_data(maskeddata == -1) = nan;
surf(nan_masked_data); shading flat;
subplot(1,3,2);
surf(log10(optimal_fit)); shading flat;
subplot(1,3,3);
residuals = maskeddata - optimal_fit;
residuals(maskeddata == -1) = nan;
surf(residuals); shading flat;
colormap(flipud(bone));



scaling_constant = 0.9;
% max(m4.disk_averages(20,20,:),)
this_blink_vals = permute(m4.disk_averages(20,20,:),[1,3,2]);
DSC_guess = [1,1];
[ pred_vals ] = trigFittingFunctions( DSC_guess, scaling_constant );

xbase = -2:0.01:2;
ybase = -2:0.01:2;
[xspace,yspace] = meshgrid(xbase,ybase);
RMSR = zeros(size(xspace));
for i = 1:length(xbase)
    i
    for j = 1:length(ybase)
        
        thisx = xspace(i,j);
        thisy = yspace(i,j);
        DSC_guess = [thisx,thisy];
        [ pred_vals ] = trigFittingFunctions( DSC_guess, scaling_constant ); 
        weight_vector = logical([1 1 1 1 1 1 1 1 1 0 1 1]);
        RMSR(i,j) = rms(this_blink_vals(weight_vector)-pred_vals(weight_vector));
    end
end

RMSR(isnan(RMSR)) = 0;

figure
contourf(xbase,ybase,RMSR);
axis equal
shading flat;

