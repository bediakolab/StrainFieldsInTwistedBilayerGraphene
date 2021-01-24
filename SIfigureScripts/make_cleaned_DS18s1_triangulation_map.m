% make_cleaned_DS18s1_triangulation_map.m
%
% Some final cleaned-up versions for publication.
%
% Nathanael Kazmierczak, 07/18/2020

fp = '/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/StrainData/AnnealedDatasetsForUpload/ChiralFixed';
fn = '05312020_chiralfixed_DS18_2_postannealing.mat';
addpath(fp);
load(fn);

m4.repairTriangulation();
m4.plotFacetedTwistAngles([0,inf],pink,false,false,true,[]);
