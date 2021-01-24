% annealing_schedule_2.m
%
% Script to help understand the best way to reconstruct the extended zone
% vector field for strain mapping.
%
% Nathanael Kazmierczak, 04/06/2020

load('04062020AnnealingStartingData.mat');
[RMSgradient_energy_storage, magnetic_energy_storage, num_vectors_changed_storage] = ...
    m4.annealDisplacementField(0.3,9,1,1,2);
m4.plotAnnealedDisplacementField();
save('04062020AnnealingSchedule2Results.mat');
clearvars m4
save('04062020AnnealingSchedule2Results_lean.mat');
