% make_explanatory_DSC_plots.m
%
% Nathanael Kazmierczak, 12/29/2019
moire_angle = 30;
[v1,v2] = getDSCBasisVectors();
trblg = TranslatedBilayerGraphene((v1+v2)/2-[0.001,0.001],moire_angle);
trblg.computeDSCField();
figh = trblg.plot();
trblg.plotDSCQuiver(figh);

% % % % 
% % % % moire_angle = 1.2;
% % % % trblg = TranslatedBilayerGraphene([0,0],moire_angle);

% % % % [v1,v2] = getDSCBasisVectors();
% % % % 
% % % % moire_angle = 1.2;
% % % % trblg = TranslatedBilayerGraphene([0,0],moire_angle);
% % % % hBN_flag = 0;
% % % % trblg.simulate(hBN_flag);
% % % % restrict_flag = 1;
% % % % 
% % % % 
% % % % trblg2 = TranslatedBilayerGraphene(v1,moire_angle);
% % % % hBN_flag = 0;
% % % % trblg2.simulate(hBN_flag);
% % % % restrict_flag = 1;
% % % % 
% % % % 
% % % % trblg3 = TranslatedBilayerGraphene((v1+v2)/2,moire_angle);
% % % % hBN_flag = 0;
% % % % trblg3.simulate(hBN_flag);
% % % % restrict_flag = 1;
% % % % 
% % % % 
% % % % trblg4 = TranslatedBilayerGraphene(v2,moire_angle);
% % % % hBN_flag = 0;
% % % % trblg4.simulate(hBN_flag);
% % % % restrict_flag = 1;
% % % % 
% % % % trblg.plotDP(restrict_flag);
% % % % trblg2.plotDP(restrict_flag);
% % % % trblg3.plotDP(restrict_flag);
% % % % trblg4.plotDP(restrict_flag);

% tblg_c = TwistedBilayerGraphene(moire_angle,0,0);
% angle = 1.2;
% distance = 45;
% tblg_1 = TwistedBilayerGraphene(moire_angle,angle,distance);
% tblg_c.plot();
% figh = tblg_c.plot();

% tblg_c.plotDSCQuiver(figh);
% tblg_c.plotDSCQuiver();
% figh2 = tblg_c.plotDSCField();
% tblg_c.assignPsuedostacking(figh2);

% tblg_c.computeDSCField();
% % tblg_c.plotDSCLattice();
% 
% tblg_1.computeDSCField();
% figh2 = tblg_1.plotDSCField();
% tblg_1.assignPsuedostacking(figh2);
% tblg_1.plotDSCLattice();


