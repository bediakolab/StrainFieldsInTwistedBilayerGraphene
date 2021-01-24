% simulate_beamspot_movement.m

xvalues = 2:0.02:8;

frac_overlap_storage = zeros(301,5);
beam_radius = 1;
reconstruction_size = 0.5;
[frac_overlap_storage(:,1)] = getBeamspotCurve(beam_radius,reconstruction_size);
beam_radius = 1;
reconstruction_size = 0.3;
[frac_overlap_storage(:,2)] = getBeamspotCurve(beam_radius,reconstruction_size);
beam_radius = 1;
reconstruction_size = 0.1;
[frac_overlap_storage(:,3)] = getBeamspotCurve(beam_radius,reconstruction_size);
beam_radius = 0.3;
reconstruction_size = 0.5;
[frac_overlap_storage(:,4)] = getBeamspotCurve(beam_radius,reconstruction_size);
beam_radius = 0.3;
reconstruction_size = 0.1;
[frac_overlap_storage(:,5)] = getBeamspotCurve(beam_radius,reconstruction_size);

figure; 
h = plot(xvalues,frac_overlap_storage,'-o','MarkerSize',3);
set(h, {'MarkerFaceColor'}, get(h,'Color')); 
legend('beam radius = 1, reconstruction size = 0.5',...
    'beam radius = 1, reconstruction size = 0.3',...
    'beam radius = 1, reconstruction size = 0.1',...
    'beam radius = 0.3, reconstruction size = 0.5',...
    'beam radius = 0.3, reconstruction size = 0.1')
title 'Real space overlap fraction'
xlabel('x')
ylabel('y')


