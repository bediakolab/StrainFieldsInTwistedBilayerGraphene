% first_script_0209_for_simulations_and_pictures.m
%
% NPK and MVW: moire simulations


moire_angle_deg = 10;
recon_angle = 0;  % total relative rotation
recon_distance = 0;  % Angstroms
use_old_reconfun = 0;  % should always be set to this.
m = TwistedBilayerGraphene(moire_angle_deg,recon_angle,recon_distance,use_old_reconfun);
% m.plot();
new_alg_flag = 1;
m.computeDSCField(new_alg_flag);
% figh = m.plotDSCField();
% m.plotDSCQuiver();
% m.plotDSCLattice();
% m.assignPsuedostacking(figh);
% m.predictBlinkingPatterns();
% m.makeBlinkingPlots(1,0);
include_hBN_flag = 0;

% [stack4D,atoms] = m.simulate(include_hBN_flag);
moire_to_match = 2;  % good number to keep (ignore)
% stacking_type = 'SP';  % pay attention to this.
stacking_type = [0.7,0.7];  % in Angstroms
tr =  TranslatedBilayerGraphene(stacking_type,moire_to_match);
include_hBN = 0;
DP = tr.simulate(include_hBN);
plotDP(DP);
