% Visualize_SP_strain.m

use_old_reconfun = 0;
recon_distance = 45;
recon_angle = 0;
moire_angle_deg = 9;
tblg = TwistedBilayerGraphene(moire_angle_deg,recon_angle,recon_distance,use_old_reconfun);
tblg.computeDSCField(0);
figh = tblg.plot();
tblg.plotDSCQuiver(figh);

recon_distance2 = 5;
recon_angle2 = 10;
tblg2 = TwistedBilayerGraphene(moire_angle_deg,recon_angle2,recon_distance2,use_old_reconfun);
tblg2.computeDSCField(0);
figh = tblg2.plot();
tblg2.plotDSCQuiver(figh);

tblg.computeNoModDSCField();
tblg2.computeNoModDSCField();
tblg.computeStrainField(0);
tblg2.computeStrainField(0);
tblg.plotStrainField();
% tblg2.plotStrainField();
% tblg.plotDSCField();
% tblg2.plotDSCField();


