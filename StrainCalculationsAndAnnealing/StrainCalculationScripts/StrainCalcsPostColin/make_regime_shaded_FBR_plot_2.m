% make_regime_shaded_FBR_plot_2.m
%
% 06/20/2020: This is the first attempt with two y-axes. Didn't work out
% all that well -- too busy. Next, try subplots.
%
% Nathanael Kazmierczak

alpha = 0.1;



addpath('/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/StrainData/StrainMappingResults/Studies/wTGVplots/MoireRemovedStandardError_06172020');
load('SPdata.mat');
figure;
subplot(2,1,1);
hold on;

% Shade according to the regime assignments
% ylimvals = get(gca,'yLim');
ylimvals = [-0.4,0.4];
xpatch = [0,0,0.6,0.6];
ypatch = [ylimvals(1),ylimvals(2),ylimvals(2),ylimvals(1)];
p = patch(xpatch,ypatch,'m');
p.FaceAlpha = alpha;
xpatch = [0.6,0.6,1.4,1.4];
ypatch = [ylimvals(1),ylimvals(2),ylimvals(2),ylimvals(1)];
p = patch(xpatch,ypatch,'c');
p.FaceAlpha = alpha;

h1 = errorbar(twist_angle_data,2*fixed_body_rotation_avs',2*fixed_body_rotation_stderrs',2*fixed_body_rotation_stderrs',xneg,xpos,'o');
% title('Fixed-body rotation in SP domains, Stderr from connected components');
xlabel('Twist angle (deg)');
ylabel('Rotation angle (deg)');
set(gca,'FontSize',14);
hold on



clearvars -except h1 alpha
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
legend([h1 h2 h3],'SP reconstruction rotation','AB reconstruction rotation','AB commensurate rotation');

% Add on ghost trendline to help the viewer see the zero mark.
xp2 = [0,1.4];
yp2 = [0,0];
plot(xp2,yp2,'-','Color',[0.8,0.8,0.8]);


% Put displacement plots.
subplot(2,1,2);
hold on;
% Shade according to the regime assignments
ylimvals = [-0.8,0];
xpatch = [0,0,0.6,0.6];
ypatch = [ylimvals(1),ylimvals(2),ylimvals(2),ylimvals(1)];
p = patch(xpatch,ypatch,'m');
p.FaceAlpha = alpha;
xpatch = [0.6,0.6,1.4,1.4];
ypatch = [ylimvals(1),ylimvals(2),ylimvals(2),ylimvals(1)];
p = patch(xpatch,ypatch,'c');
p.FaceAlpha = alpha;

[ v1, v2, hexagon_lattice_constant ] = getDSCBasisVectors();
a = hexagon_lattice_constant;
thetam = deg2rad(twist_angle_data);
r = a./(thetam*2*sqrt(3));
thetaAB = deg2rad(2*fixed_body_rotation_avs');  % multiply by two so it is top and bottom on the domain.
u = thetaAB.*r;
h4 = scatter(twist_angle_data,u);
hold on
xlabel('Twist angle (degrees)');
ylabel('Displacement (Angstrom)');

% Put limiting displacement value 
xpd = [1.4,0.05];
ypd = [-a/sqrt(3)/2,-a/sqrt(3)/2];
h5 = plot(xpd,ypd,'-');
legend([h4,h5],'AB-induced SP displacement','Perfect soliton SP displacement');

