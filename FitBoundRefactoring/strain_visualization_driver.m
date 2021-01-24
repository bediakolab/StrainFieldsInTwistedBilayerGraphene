% strain_visualization_driver.m

tblg = TwistedBilayerGraphene(1.2,0.5,40,0);
tblg.computeNoModDSCField();
tblg.plotNoModDSCField();
tblg.computeStrainField(0);
tblg.plotStrainField();

