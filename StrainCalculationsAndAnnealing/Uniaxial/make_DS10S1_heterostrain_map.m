% make_DS10S1_heterostrain_map.m
%
% Nathanael Kazmierczak, 10/23/2020

addpath('/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/StrainData/PreparedDatasets/ChiralFixed');
load('05312020_chiralfixed_DS10_forannealing.mat');

m4.repairTriangulation();

angle_guess = 0.63;  % right on for DS8S1
average_angle_range = [];
[average_rms_residual,nconvergedstor] = m4.fitUniaxialStrain(average_angle_range,angle_guess);
% make_plots = 1;
% moire_binedges = [];
% strainangle_binedges = [];
% strainval_binedges = 0:0.025:0.45;
% [av_moire_angle,std_moire_angle,min_moire_angle,max_moire_angle,...
%                 av_strain_percent,std_strain_percent,min_strain_percent,max_strain_percent,copypaste] = ...
%                 m4.getUniaxialStrainStatistics(make_plots,moire_binedges,strainangle_binedges,strainval_binedges);

removeLines = true;     
cmap = magma;
m4.plotFacetedUniaxialStrain(cmap,removeLines);

close(gcf)

xlabel('x (nm)');
ylabel('y (nm)');
