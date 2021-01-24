% triangulation_script.m
%
% To layer histograms for Figure 2b.
% Nathanael Kazmierczak, 07/05/2020

addpath('/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/GeometryFitting/TriangulationPlots');
load('Triangulation_data_for_plotting_07052020.mat');
load('DS2S1_triangulation_data.mat');

figure;
h1 = histogram(DS4S2_twistdata,1.06:0.02:1.3,'Normalization','probability');
hold on;
h2 = histogram(DS5S1_twistdata,1.14:0.02:1.34,'Normalization','probability');
h3 = histogram(DS8S1_twistdata,1.26:0.02:1.48,'Normalization','probability');
xlabel('Twist angle (deg)');
ylabel('Area fraction');
legend('\Theta_m = 1.19^o','\Theta_m = 1.23^o','\Theta_m = 1.37^o');
title('Version 1');

figure;
h1 = histogram(DS2S1_twistdata,0.94:0.02:1.14,'Normalization','probability');
hold on;
h2 = histogram(DS5S1_twistdata,1.14:0.02:1.34,'Normalization','probability');
h3 = histogram(DS8S1_twistdata,1.26:0.02:1.48,'Normalization','probability');
xlabel('Twist angle (deg)');
ylabel('Area fraction');
legend('\Theta_m = 1.03^o','\Theta_m = 1.23^o','\Theta_m = 1.37^o');
title('Version 2');

figure;
h1 = histogram(DS2S1_twistdata,0.94:0.02:1.14,'Normalization','probability');
hold on;
h2 = histogram(DS4S2_twistdata,1.06:0.02:1.3,'Normalization','probability');
h3 = histogram(DS8S1_twistdata,1.26:0.02:1.48,'Normalization','probability');
xlabel('Twist angle (deg)');
ylabel('Area fraction');
legend('\Theta_m = 1.03^o','\Theta_m = 1.19^o','\Theta_m = 1.37^o');
title('Version 3');
