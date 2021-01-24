% plotAAshearstrain.m
%
% Nathanael Kazmierczak, 05/28/2020

error('This is an old script that used incorrect conventions for computing strain');

gxy_means = [1.5834
1.2324
1.2916
1.0036
2.5335
4.4898
4.5399
4.2811];

moire_angles = [1.0325
1.2276
1.2343
1.3672
0.6272
0.1203
0.1355
0.2609];

gxy_stds = [0.2336
0.2023
0.2486
0.2068
0.2742
0.55
0.4793
0.6247];

theta_stds = [0.0333,0.024,0.0284,0.0315,0.0189,0.00071638,0.0067,0.0248]';

figure;
yneg = gxy_stds;
ypos = gxy_stds;
xneg = theta_stds;
xpos = theta_stds;
errorbar(moire_angles,gxy_means,yneg,ypos,xneg,xpos,'o');
title('gxy total shear strain in AA domains');
xlabel('Twist angle (degrees)');
ylabel('Strain %');

