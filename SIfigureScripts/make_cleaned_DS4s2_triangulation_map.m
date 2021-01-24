% make_cleaned_DS4s2_triangulation_map.m
%
% Some final cleaned-up versions for publication.
%
% Nathanael Kazmierczak, 07/18/2020

fp = '/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/StrainData/PreparedDatasets/Session2';
fn = '06222020_DS4S2_prepared_for_annealing.mat';
addpath(fp);
load(fn);

m4.plotFacetedTwistAngles([0,inf],pink,false,false,true,[]);
