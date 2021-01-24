% make_cleaned_DS12s2_triangulation_map.m
%
% Some final cleaned-up versions for publication.
%
% Nathanael Kazmierczak, 07/18/2020

fp = '/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/StrainData/AnnealedDatasetsForUpload/Session2';
fn = '06222020_DS12S2_annealed.mat';
addpath(fp);
load(fn);

m4.plotFacetedTwistAngles([0,inf],pink,false,false,true,[]);
