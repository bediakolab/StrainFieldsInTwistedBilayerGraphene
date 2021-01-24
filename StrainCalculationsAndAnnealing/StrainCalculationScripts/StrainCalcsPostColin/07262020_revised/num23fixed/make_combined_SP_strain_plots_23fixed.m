% make_combined_SP_strain_plots_23fixed.m
%
% Script for combining the results of three SP calculations with different
% tensor rotations.
%
% Nathanael Kazmierczak, 07/27/2020

filenames = {'SP1_23fixed.mat';
             'SP2_23fixed.mat';
             'SP3_23fixed.mat'};

load(filenames{1});
sxy_avs_stor(:,:,1) = sxy_avs;
sxy_stderrs_stor(:,:,1) = sxy_stderrs;
syx_avs_stor(:,:,1) = syx_avs;
syx_stderrs_stor(:,:,1) = syx_stderrs;
gammaxy_avs_stor(:,:,1) = principal_shear_avs;
gammaxy_stderrs_stor(:,:,1) = principal_shear_stderrs;
FBR_stor(:,:,1) = fixed_body_rotation_avs;
FBR_stderrs(:,:,1) = fixed_body_rotation_stderrs;

load(filenames{2});
sxy_avs_stor(:,:,2) = sxy_avs;
sxy_stderrs_stor(:,:,2) = sxy_stderrs;
syx_avs_stor(:,:,2) = syx_avs;
syx_stderrs_stor(:,:,2) = syx_stderrs;
gammaxy_avs_stor(:,:,2) = principal_shear_avs;
gammaxy_stderrs_stor(:,:,2) = principal_shear_stderrs;
FBR_stor(:,:,2) = fixed_body_rotation_avs;
FBR_stderrs(:,:,2) = fixed_body_rotation_stderrs;

load(filenames{3});
sxy_avs_stor(:,:,3) = sxy_avs;
sxy_stderrs_stor(:,:,3) = sxy_stderrs;
syx_avs_stor(:,:,3) = syx_avs;
syx_stderrs_stor(:,:,3) = syx_stderrs;
gammaxy_avs_stor(:,:,3) = principal_shear_avs;
gammaxy_stderrs_stor(:,:,3) = principal_shear_stderrs;
FBR_stor(:,:,3) = fixed_body_rotation_avs;
FBR_stderrs(:,:,3) = fixed_body_rotation_stderrs;


FBRplot = mean(FBR_stor,3);
sxyplot = mean(sxy_avs_stor,3);
syxplot = mean(syx_avs_stor,3);
gammaxyplot = mean(gammaxy_avs_stor,3);
FBRstderrplot = mean(FBR_stderrs,3);
sxy_stderrplot = mean(sxy_stderrs_stor,3);
syx_stderrplot = mean(syx_stderrs_stor,3);
gammaxy_stderrplot = mean(gammaxy_stderrs_stor,3);

%% Make plots
figure
errorbar(twist_angle_data_use,2*FBRplot',2*2*FBRstderrplot',2*2*FBRstderrplot',twist_angle_stderrs_use,twist_angle_stderrs_use,'-o');
title('SP FBR averaged over all three SPs');
xlabel('Twist angle (degrees)');
ylabel('Rotation angle (degrees)');
set(gca,'FontSize',14);
legend('SP rotation (2 layers), 95% CI error bars');

figure
errorbar(twist_angle_data_use,sxyplot',2*sxy_stderrplot',2*sxy_stderrplot',twist_angle_stderrs_use,twist_angle_stderrs_use,'-o');
hold on
errorbar(twist_angle_data_use,syxplot',2*syx_stderrplot',2*syx_stderrplot',twist_angle_stderrs_use,twist_angle_stderrs_use,'-o');
hold on
errorbar(twist_angle_data_use,gammaxyplot',2*gammaxy_stderrplot',2*gammaxy_stderrplot',twist_angle_stderrs_use,twist_angle_stderrs_use,'-o');
title('Shear strain in SP domains, Stderr from connected components');
xlabel('Twist angle (degrees)');
ylabel('Strain %')
set(gca,'FontSize',14);
legend('sxy','syx','Gamma max');


