% gaussianFBRprofile.m

addpath('/Users/nathanaelkazmierczak/Dropbox/NPKHardDriveBackupFall2017/workspace/BediakoResearch/NPK_MatlabDataGeneral/Reports/LocalFigureDrafts/ExtendedData/07212020/Donotupload');
openfig('DS26FBRprofile_init_AB_AB.fig');
f1 = gcf;
axobjs = f1.Children;
dataObjs = axobjs.Children;
xvals =  dataObjs.XData;
yvals =  dataObjs.YData;
newy = yvals * 2;  % Conversion to interlayer as desired. 
[~,maxind] = max(newy);
xmax = xvals(maxind);
newx = xvals - xmax;

% similar operation for the SP one.
openfig('DS26FBRprofile_init_SP_SP.fig');
f2 = gcf;
axobjs = f2.Children;
dataObjs = axobjs.Children;
xvals =  dataObjs.XData;
yvals =  dataObjs.YData;
newy2 = yvals * 2;  % Conversion to interlayer as desired. 
[~,maxind] = max(newy2);
xmax = xvals(maxind);
newx2 = xvals - xmax;

% manual shift to align the x axes
newx2 = newx2 + 0.75;

figure;
plot(newx,newy,'-');
hold on
plot(newx2,newy2);
xlabel('Distance from AA center (nm)');
ylabel('Interlayer fixed-body rotation (degrees)');
legend('AB-AA-AB','SP-AA-SP');

