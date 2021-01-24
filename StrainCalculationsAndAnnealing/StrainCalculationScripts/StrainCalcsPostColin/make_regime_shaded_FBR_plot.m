% make_regime_shaded_FBR_plot.m
%
% 06/20/2020: This is the first attempt with two y-axes. Didn't work out
% all that well -- too busy. Next, try subplots.
%
% Nathanael Kazmierczak


addpath('/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/StrainData/StrainMappingResults/Studies/wTGVplots/MoireRemovedStandardError_06172020');
load('SPdata.mat');
figure;
yyaxis left
h1 = errorbar(twist_angle_data,2*fixed_body_rotation_avs',2*fixed_body_rotation_stderrs',2*fixed_body_rotation_stderrs',xneg,xpos,'o','Color',[0,0,0.6]);
title('Fixed-body rotation in SP domains, Stderr from connected components');
xlabel('Twist angle (degrees)');
ylabel('Rotation angle (degrees)');
set(gca,'FontSize',14);
hold on

clearvars -except h1
load('ABdata.mat');
yneg = fixed_body_rotation_stderrs';
ypos = fixed_body_rotation_stderrs';
xneg = twist_angle_stds;
xpos = twist_angle_stds;
h2 = errorbar(twist_angle_data,2*fixed_body_rotation_avs',2*yneg,2*ypos,xneg,xpos,'o');

% Add on linear trendline for small theta m
xp = [0,0.3];
yp = [0,-0.3];
h3 = plot(xp,yp,'-');


% Put displacement plots.
yyaxis right
[ v1, v2, hexagon_lattice_constant ] = getDSCBasisVectors();
a = hexagon_lattice_constant;
thetam = deg2rad(twist_angle_data);
r = a./(thetam*2*sqrt(3));
thetaAB = deg2rad(2*fixed_body_rotation_avs');  % multiply by two so it is top and bottom on the domain.
u = thetaAB.*r;
h4 = scatter(twist_angle_data,u);
ylabel('Saddle Point Displacement (nm)');

% Put limiting displacement value 
xpd = [1.4,0.05];
ypd = [-a/sqrt(3)/2,-a/sqrt(3)/2];
h5 = plot(xpd,ypd,'-');


legend([h1 h2 h3 h4 h5],'SP reconstruction rotation','AB reconstruction rotation','AB commensurate rotation','Induced SP displacement','Perfect soliton SP displacement');


