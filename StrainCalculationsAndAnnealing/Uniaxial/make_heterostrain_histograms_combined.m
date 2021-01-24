% make_heterostrain_histograms_combined.m
%
% Nathanael Kazmierczak, 10/23/2020

addpath('/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/StrainData/AnnealedDatasetsForUpload/ChiralFixed');
load('05312020_chiralfixed_DS8_postannealing.mat');
m4.repairTriangulation(false);
angle_guess = 1.37;  % right on for DS8S1
average_angle_range = [];
[average_rms_residual,nconvergedstor] = m4.fitUniaxialStrain(average_angle_range,angle_guess);
DS8S1_heterostrain = m4.uniaxial_strain_percents;
clear m4

load('05312020_chiralfixed_DS4_forannealing.mat');
m4.repairTriangulation(false);
angle_guess = 1.23;  % DS4S1
average_angle_range = [];
[average_rms_residual,nconvergedstor] = m4.fitUniaxialStrain(average_angle_range,angle_guess);
DS4S1_heterostrain = m4.uniaxial_strain_percents;
clear m4

load('05312020_chiralfixed_DS2_forannealing.mat');
m4.repairTriangulation(false);
angle_guess = 1.03;  % DS4S1
average_angle_range = [];
[average_rms_residual,nconvergedstor] = m4.fitUniaxialStrain(average_angle_range,angle_guess);
DS2S1_heterostrain = m4.uniaxial_strain_percents;
clear m4

addpath('/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/StrainData/PreparedDatasets/Session2');
load('06222020_DS4S2_prepared_for_annealing.mat');
m4.repairTriangulation(false);
angle_guess = 1.19;  % DS8S2
average_angle_range = [];
[average_rms_residual,nconvergedstor] = m4.fitUniaxialStrain(average_angle_range,angle_guess);
DS4S2_heterostrain = m4.uniaxial_strain_percents;
clear m4


bins = 0:0.05:0.6;
figure;
histogram(DS8S1_heterostrain,bins);
hold on;
histogram(DS4S1_heterostrain,bins);
histogram(DS4S2_heterostrain,bins);
histogram(DS2S1_heterostrain,bins);
xlabel('Heterostrain %');
ylabel('Counts');

legend('DS8S1 \Theta=1.37^{\circ}','DS4S1 \Theta=1.23^{\circ}','DS4S2 \Theta=1.19^{\circ}','DS2S1 \Theta=1.03^{\circ}');

