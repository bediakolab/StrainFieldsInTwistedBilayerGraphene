% testCutROI.m
%
% Script for testing the function cutROI.m

addpath('/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/matlab4DSTEM');
load('08272019testingCutROI.mat')

figure; image(uint8(firstDP),'CDataMapping','Scaled'); colormap pink;
index_choice = [0,1];
window_width = 15;
% first_graphene_disks = fliplr(first_graphene_disks);  % because they were saved as x,y rather than I,J
[ roicut, Irangecut, Jrangecut ] = cutROI( firstDP, first_graphene_disks, hk_graphene1, index_choice, window_width );
Irangecut
Jrangecut
figure; image(uint8(roicut),'CDataMapping','Scaled'); colormap pink;

