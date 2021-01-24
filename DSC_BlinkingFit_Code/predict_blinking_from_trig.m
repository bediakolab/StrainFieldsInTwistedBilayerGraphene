
% predict_blinking_from_trig.m
%
% Driver script for making blinking predictions for several types of
% reconstruction.
%
% Nathanael Kazmierczak, 12/28/2019

addpath(genpath(pwd));

moire_angle = 1.2;
recon_angle = 0.8;
recon_distance = 45;

tblg1 = TwistedBilayerGraphene(moire_angle,recon_angle,recon_distance);
tblg_c = TwistedBilayerGraphene(moire_angle,0,0);
tic;
tblg1.computeDSCField();
tblg_c.computeDSCField();
toc

% tblg.plot();

figh = tblg1.plotDSCField();
figh = tblg1.assignPsuedostacking(figh);
fighc = tblg_c.plotDSCField();
fighc = tblg_c.assignPsuedostacking(fighc);
% % % tblg.plotDSCLattice();
% % % figh = tblg.assignPsuedostacking(figh);
% % recon_angle = 0.8;
% % recon_distance = 40;
% % tblg2 = TwistedBilayerGraphene(moire_angle,recon_angle,recon_distance);
% % recon_angle = 5;
% % recon_distance = 40;
% % tblg3 = TwistedBilayerGraphene(moire_angle,recon_angle,recon_distance);
% % 
% % tblg1.computeDSCField();
% % tblg2.computeDSCField();
% % tblg3.computeDSCField();
% % 
subplot_flag = 1;
include_redundant_flag = 0;
tblg1.makeBlinkingPlots(subplot_flag,include_redundant_flag);
tblg_c.makeBlinkingPlots(subplot_flag,include_redundant_flag);
% % tblg2.makeBlinkingPlots(subplot_flag,include_redundant_flag);
% % tblg3.makeBlinkingPlots(subplot_flag,include_redundant_flag);
