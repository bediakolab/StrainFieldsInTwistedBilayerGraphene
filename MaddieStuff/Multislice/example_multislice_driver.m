% example_multislice_driver.m
%
% This script sketches out how NPK used the TwistedBilayerGraphene.m class
% in support of developing an interferometry fitting function. The idea is
% to simulate diffraction patterns for various atomic structures
% characterized by displacement vectors, extract the intensities of regions
% of interest, plot as a function of displacement vector, and use that to
% determine an appropriate fitting model.
%
% Upon successful execution of this script, three figures will appear.
% Figure 1 displays the atomic coordinates of the graphene bilayer (not
% twisted!). Figure 2 displays the diffraction pattern with an additional
% hBN encapsulating layer, while Figure 3 displays the diffraction pattern
% with only the graphene and central beam considered. 
%
% Nathanael Kazmierczak, 08/13/2020


moire_cell_to_match = 2;  % 10 degrees gives a small unit cell, makes for a fast simulation.
DSC_vector_or_type = [0.4,0.75]; % cartesian displacement coordinates
trblg = TranslatedBilayerGraphene(DSC_vector_or_type,moire_cell_to_match);
trblg.plot();

include_hBN_flag = true;
extra_spacing_factor = 0;  
% I tested the spacing once to see if z-axis information could be obtained
% from intensity. The answer seemed to be no -- Hamish's "projected
% potential" concept in the derivation of the fitting function (see Methods
% of TBLG interferometry paper) is probably why.
[stack4D,atoms] = trblg.simulate(include_hBN_flag,extra_spacing_factor);
restrict_flag = false;
trblg.plotDP(restrict_flag);
caxis([0,1e-5]);  % make diffraction spots visible around center beam.

include_hBN_flag = false;
[stack4D,atoms] = trblg.simulate(include_hBN_flag,extra_spacing_factor);
trblg.plotDP(restrict_flag);
caxis([0,1e-5]);
