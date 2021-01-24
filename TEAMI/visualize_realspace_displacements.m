% visualize_realspace_displacements.m

moire_angle_deg = 1.1;
recon_angle = 3;
recon_distance = 50;
use_old_reconfun = 0;
tblg = TwistedBilayerGraphene(moire_angle_deg,recon_angle,recon_distance,use_old_reconfun);
newalgflag = 0;
tblg.computeDSCField(newalgflag);
tblg.computeNoModDSCField();
figh = tblg.plotDSCField();
% tblg.plotNoModDSCField();

tblg.makeInteractiveDisplacementCircles();


