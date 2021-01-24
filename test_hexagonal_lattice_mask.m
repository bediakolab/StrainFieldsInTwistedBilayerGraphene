% test_hexagonal_lattice_mask.m
%
% Script for testing the new hexagonal lattice masking method.
% Function is makeHexagonalLatticeMask.m
% Nathanael Kazmierczak, Nov 2019

% Suck in an example diffraction pattern that we are using for the
% ptychography.


addpath('/Volumes/NPKResearch/BediakoLab/Maddie4DSTEM/New2019_1023');
addpath('/Volumes/NPKResearch/BediakoLab/Maddie4DSTEM/New2019_1023/Dataset19');
addpath('/Volumes/NPKResearch');
addpath('/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/matlab4DSTEM');
filename = '19_MV_1p5_9_30x30_ss1nm_1s_spot9_alpha=1p75_bin2_cl=130_60kV_defocus=m400NPK_cropped.h5';
dataset_name = '/4DSTEM_experiment/data/datacubes/datacube_0/data';
if ~exist('data','var')
    % Size is [1024 1024 19 14]
    info = h5info(filename);
    starting = [1,1,1,1];
    ending = [1024,1024,1,1];
    DP1 = double(h5read(filename,dataset_name,starting,ending));
    disp('Data has been loaded from USB drive.');
else
    disp('Using data already loaded into workspace!');
end

figure;
exponent = 0.5;
imagesc(DP1.^exponent);
axis square;
colormap jet
colorbar

radius = 45;
[ masked_image, mask, outlattice, beamstop_mask, outbasis ] = makeHexagonalLatticeMask( DP1, radius );

radius = 35;
[ masked_image, mask, outlattice, beamstop_mask ] = ...
    makeHexagonalLatticeMask( DP1, radius, outlattice, beamstop_mask );

radius = 40;
[ masked_image, mask, outlattice, beamstop_mask ] = ...
    makeHexagonalLatticeMask( DP1, radius, outlattice, beamstop_mask );


