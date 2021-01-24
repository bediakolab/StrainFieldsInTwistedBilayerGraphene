% make_atom_displacements_plot.m
%
% Nathanael Kazmierczak, 10/30/2020

if ~exist('m4','var')
    addpath('/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/StrainData/AnnealedDatasetsForUpload/ChiralFixed');
    load('05312020_chiralfixed_DS4_postannealing.mat');
end
dfield = m4.annealed_dfield;
spacing = m4.scan_stepsize;

% moire_angle_deg = 3;
moire_angle_deg = 0;
recon_struct.type = 'none';
heterostrain_struct = [];
custom_cell_size = [1000,1000];
% custom_cell_size = [50,50];
tblg = TwistedBilayerGrapheneAugmented(moire_angle_deg,recon_struct,heterostrain_struct,custom_cell_size);
% tblg.invertCoords();
tblg.splineDisplacement(dfield,spacing);
tblg.plot();
% tblg.plotVolumetric();
