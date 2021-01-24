% team1_hotspots.m
%
% Script for interactively probing the reasons for the hotspots.
%
% 03/06/2020
% Nathanael Kazmierczak

if false
    load('TEAMI_Day1_DS4_objectdata.mat');
    m4.makeBlinkingPlots([],0,'flat');
    m4.makeDisplacementMapsFromBlinking(0,'flat',fire,0);
    m4.makeInteractiveDisplacementDarkFieldPlot();
    m4.makeInteractiveInverseDarkFieldPlot();
end

if true
    load('TEAMI_Day1_DS7_objectdata.mat');
    m4.makeBlinkingPlots([],0,'flat');
    m4.makeDisplacementMapsFromBlinking(0,'flat',fire,0);
    m4.makeInteractiveDisplacementDarkFieldPlot();
    m4.makeInteractiveInverseDarkFieldPlot();
end

