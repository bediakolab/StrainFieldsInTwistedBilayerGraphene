% make_cleaned_DS3s2_displacement_map.m
%
% Some final cleaned-up versions for publication.
%
% Nathanael Kazmierczak, 07/18/2020

fp = '/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/StrainData/AnnealedDatasetsForUpload/Session2';
fn = '06222020_DS3S2_annealed.mat';
addpath(fp);
load(fn);

stringcell = {'soft multistart','amplitude',[0.3,5];'median outlier','amplitude',[5,0.3]};
filterstruct = buildFilterStruct(stringcell);
m4.makeCustomDisplacementColorPlot([],[],filterstruct,[],1,0,0,0);

