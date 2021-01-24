% positionfit_0203_DS8_testing.m
%
% Trying both the blinking and the position strain analysis methods.

fit_blinking_flag = false;
fit_position_flag = true;

filename = '8_80x80_ss=1nm_C2=40um_alpha=1mrad_spot7_100ms_CL=130 bin=4_rot=90_60kV.h5';
probe_filename = []; % Replace this once Maddie uploads it. '60kV_1mrad_spot9_0p5s_130CL_defocus = 0.dm3';
currentdir = pwd;
cd ..
addpath(genpath(pwd));
cd(currentdir);
addpath('/Volumes/NPKResearch/BediakoLab/Maddie4DSTEM/20200203');
results_path = '/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/PositionFit/Results/20200203_dataset8';

% % disks_name = '9_40x40_ss=0p5nm_C2=40um_alpha=1mrad_spot9_1s_CL=130_bin=2_scan35_60kV_NPK_Normalized__integration_disks.mat';
% disks_name = 'ds9disks.mat';
% disks_name = '7_80x80_ss=2nm_C2=40um_alpha=1mrad_spot7_200ms_CL=130 bin=4_60kV_NPKnormalized__integration_disks.mat';
disks_name = [];
stepsize = 1; % given the RS binning of 2
m4 = FourDSTEM_Analysis_Engine(filename,disks_name,results_path,stepsize,probe_filename);
% explicit_chunks = [15,16]

if fit_position_flag
    subtract_power_law = 0;
    [average_BN_probe,xcoords,ycoords] = m4.setKernelFromHBN(subtract_power_law);
    m4.plotBNprobe();
    m4.setPositionFitLattices();
    use_BN_kernel = 1;
    m4.fitDiskPositions(use_BN_kernel);
end

if fit_blinking_flag
    m4.fitElasticPowerLaw([]);
    m4.plotPowerLawMask();
    subtract_power_law_flag = 1;
    m4.setBlinkingDisks();
    m4.integrateDisks(subtract_power_law_flag);
    saveplotflag = 1;
    trig_prefactor = 0.95;
    individual_residuals = 1;
    shadingtype = 'flat';
    m4.makeBlinkingPlots([],saveplotflag,shadingtype);
    
    
    m4.fitBlinking(trig_prefactor); % no override
    median_filter_thresholds = [1,0.5];
    m4.makeDisplacementMapsFromBlinking(saveplotflag,median_filter_thresholds,individual_residuals);
    colormaptype = 'gray';
    trimwidth = 0;
    m4.makeStrainMapsFromBlinking(saveplotflag,shadingtype,colormaptype,trimwidth,median_filter_thresholds);
    m4.saveObject();
    m4.makeQuiverPlotFromBlinking(saveplotflag);
    m4.makeInteractiveDisplacementDarkFieldPlot();
end
