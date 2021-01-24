% make_cleaned_DS26s1_displacement_map.m
%
% Some final cleaned-up versions for publication.
%
% Nathanael Kazmierczak, 07/18/2020

fp = '/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/StrainData/AnnealedDatasetsForUpload/ChiralFixed';
fn = '05312020_chiralfixed_noprefactor_DS26_postannealing.mat';
addpath(fp);
load(fn);

stringcell = {'soft multistart','amplitude',[0.3,9];'median outlier','amplitude',[5,0.7]};
filterstruct = buildFilterStruct(stringcell);
m4.makeCustomDisplacementColorPlot([],[],filterstruct,[],1,0,0,0);

