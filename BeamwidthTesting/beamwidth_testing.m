% beamwidth_testing.m
%
% Simulation for determining if the displacement vector holes can be
% attributed to the finite beam width.


moire_angle_deg = 0.5;
recon_angle = 0.5;
recon_distance = 50;
use_old_reconfun = 0;
tblg = TwistedBilayerGraphene(moire_angle_deg,recon_angle,recon_distance,use_old_reconfun);
new_alg_flag = 1;
tblg.computeDSCField(new_alg_flag);
% figh = tblg.plotDSCField();
% tblg.assignPsuedostacking(figh);
% tblg.plotDSCLattice();
stepsize = 5;
averaging_radius = 50;
tblg.plotAveragedDSCLattice(stepsize,averaging_radius);


