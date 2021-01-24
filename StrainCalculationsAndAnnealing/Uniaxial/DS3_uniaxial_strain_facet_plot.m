% DS3_uniaxial_strain_facet_plot.m
%
% For generating the uniaxial model strain % plot with triangulated facets.
%
% Nathanael Kazmierczak, 07/05/2020

addpath('/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/StrainData/AnnealedDatasetsForUpload/ChiralFixed');
load('06012020_chiralfixed_tear_DS3_postannealing.mat');
% m4.repairTriangulation();
average_angle_range = [0,inf];
angle_guess = 1.0;
[average_rms_residual,nconvergedstor] = m4.fitUniaxialStrain(average_angle_range,angle_guess);

% [av_moire_angle,std_moire_angle,min_moire_angle,max_moire_angle,...
%                 av_strain_percent,std_strain_percent,min_strain_percent,max_strain_percent,copypaste] = ...
%                 getUniaxialStrainStatistics(obj,make_plots,moire_binedges,strainangle_binedges,strainval_binedges);

removeLines = true;
cmap = viridis;
m4.plotFacetedUniaxialStrain(cmap,removeLines);
close  % Get rid of heterostrain angle plot, which does not currently work.
title viridis
cmap = plasma;
m4.plotFacetedUniaxialStrain(cmap,removeLines);
close
title plasma
cmap = magma;
m4.plotFacetedUniaxialStrain(cmap,removeLines);
close
title magma
cmap = inferno;
m4.plotFacetedUniaxialStrain(cmap,removeLines);
close
title inferno

            