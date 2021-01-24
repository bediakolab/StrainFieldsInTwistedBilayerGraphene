% plotAAfixed_body_rotation.m
%
% Doing this the correct way after reading more about continuum mechanics.
% 
% Nathanael Kazmierczak, 05/31/2020

fixed_body_rotation_means = rad2deg([1.5834
1.2324
1.2916
1.0036
2.5335
4.4898
4.5399
4.2811]/100);

moire_angles = [1.0325
1.2276
1.2343
1.3672
0.6272
0.1203
0.1355
0.2609];

fixed_body_rotation_stds = rad2deg([0.2336
0.2023
0.2486
0.2068
0.2742
0.55
0.4793
0.6247]/100);

theta_stds = [0.0333,0.024,0.0284,0.0315,0.0189,0.00071638,0.0067,0.0248]';

total_twist_angle = moire_angles + fixed_body_rotation_means;

figure;
yneg = fixed_body_rotation_stds;
ypos = fixed_body_rotation_stds;
xneg = theta_stds;
xpos = theta_stds;
errorbar(moire_angles,fixed_body_rotation_means,yneg,ypos,xneg,xpos,'o');
hold on
plot(moire_angles,moire_angles,'-');
errorbar(moire_angles,total_twist_angle,yneg,ypos,xneg,xpos,'o');
title('Fixed body rotation in AA domains');
xlabel('Moire twist angle (degrees)');
ylabel('Twist angle (degrees)');
legend('Reconstruction rotation (net, over two layers)','Moire rotation','Total rotation');

