% initial_driver.m

% addpath('/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/Code_from_Tadmor/matlab_GRMED/NPKnotes');
% addpath('/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/matlab4DSTEM');
mydir = pwd;
cd ..
addpath(genpath(pwd))
cd(mydir);

moire_angle = 1.2;
recon_angle = 0.8;
recon_distance = 40;
% recon_angle = 0;
% recon_distance = 0;
tblg = TwistedBilayerGraphene(moire_angle,recon_angle,recon_distance);
tblg.plot();
% figure;
% hold on;
tblg.computeDSCField();
figh = tblg.plotDSCField();
tblg.plotDSCLattice();
figh = tblg.assignPsuedostacking(figh);

% tblg02 = TwistedBilayerGraphene(moire_angle,0,0);
% tblg02.plot();
% figure;
% hold on;
% tblg02.computeDSClattice();
% tblg02.plotDSC();

% figure; contourf(DSCamp,50,'LineStyle','None');
