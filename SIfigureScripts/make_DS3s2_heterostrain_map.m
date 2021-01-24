% make_DS3s2_heterostrain_map.m
%
% Some final cleaned-up versions for publication.
%
% Nathanael Kazmierczak, 07/18/2020

fp = '/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/StrainData/AnnealedDatasetsForUpload/Session2';
fn = '06222020_DS3S2_annealed.mat';
addpath(fp);
load(fn);


m4.fitUniaxialStrain([0,inf],1.2);
[av_moire_angle,std_moire_angle,min_moire_angle,max_moire_angle,...
                av_strain_percent,std_strain_percent,min_strain_percent,max_strain_percent,copypaste] = ...
                m4.getUniaxialStrainStatistics(false,[],[],[])
m4.plotFacetedUniaxialStrain(inferno,true);
