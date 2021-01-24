% fitting_schematic_DP.m

fp = '/Volumes/NPKResearch/BediakoLab/Maddie4DSTEM/20200226/Dataset4';
fn = 'TEAM1Day1Dataset4DiffractionSI.h5';

fp2 = '/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/TEAM1Data/withPrefactorFit/DS4/AAstatistics';
fn2 = 'DS4_prefactorfit_objectdata.mat';

addpath(fp);
addpath(fp2);
if ~exist('m4','var');
    load(fn2);
end
m4.filename = fn;
% Now it can access files off of my USB drive.
idx1 = 60;
idx2 = 103;
% idx1 = 61;
% idx2 = 104;

make_histograms_flag = false;
remove_hBN = false;
remove_Graphene = false;
remove_Beamstop = false;
correct_DP = false;
m4.plotSingleDP(idx1,idx2,make_histograms_flag,remove_hBN,remove_Graphene,remove_Beamstop,correct_DP);
close
caxis([0,16]);
colormap(gray);

remove_hBN = false;
remove_Graphene = false;
remove_Beamstop = true;
m4.plotSingleDP(idx1,idx2,make_histograms_flag,remove_hBN,remove_Graphene,remove_Beamstop,correct_DP);
close
caxis([0,16]);
colormap(gray);

remove_hBN = true;
remove_Graphene = false;
remove_Beamstop = true;
m4.plotSingleDP(idx1,idx2,make_histograms_flag,remove_hBN,remove_Graphene,remove_Beamstop,correct_DP);
close
caxis([0,16]);
colormap(gray);

exponent = 1;
removeHBN = false;
removeGraphene = false;
removeBeamstop = false;
subtractScattering = false;
subtractRadial = false;
subtractCartesianInterpolation = false;
subtractPolarInterpolation = false;
m4.plotBaselinedData(exponent,removeHBN,removeGraphene,removeBeamstop,subtractScattering,subtractRadial,subtractCartesianInterpolation,subtractPolarInterpolation);
colormap(fire);
caxis([0,4]);

exponent = 1;
removeHBN = true;
removeGraphene = false;
removeBeamstop = false;
subtractScattering = true;
subtractRadial = true;
subtractCartesianInterpolation = false;
subtractPolarInterpolation = true;
m4.plotBaselinedData(exponent,removeHBN,removeGraphene,removeBeamstop,subtractScattering,subtractRadial,subtractCartesianInterpolation,subtractPolarInterpolation);
colormap(fire);
caxis([0,2.5]);

exponent = 1;
removeHBN = true;
removeGraphene = false;
removeBeamstop = false;
subtractScattering = false;
subtractRadial = false;
subtractCartesianInterpolation = false;
subtractPolarInterpolation = false;
m4.plotBaselinedData(exponent,removeHBN,removeGraphene,removeBeamstop,subtractScattering,subtractRadial,subtractCartesianInterpolation,subtractPolarInterpolation);
colormap(fire);
caxis([0,4]);
