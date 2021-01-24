% find_DSC_blinking_functional_form2.m

if ~exist('averages','var')
    load('12272019_trblg_linesearches2_data2.mat','averages');
    load('12272019_trblg_linesearches2_data2.mat','lines');
end

DSC_values = lines{1};
theseaverages = averages(:,:,1);
disk_num = [1:3,7:9];
figh = figure;
makeProjectedBlinkingPlot( disk_num, theseaverages, DSC_values);
% 
% disk_num = 2;
% makeProjectedBlinkingPlot( disk_num, theseaverages, DSC_values, figh);
% 
% disk_num = 3;
% makeProjectedBlinkingPlot( disk_num, theseaverages, DSC_values, figh);
% 
% disk_num = 7;
% makeProjectedBlinkingPlot( disk_num, theseaverages, DSC_values, figh);
% 
% disk_num = 8;
% makeProjectedBlinkingPlot( disk_num, theseaverages, DSC_values, figh);
% 
% disk_num = 9;
% makeProjectedBlinkingPlot( disk_num, theseaverages, DSC_values, figh);

% 
% [ dict_function, v ] = getLatticePlaneNormalVectors();
% 
% line1 = lines{1};
% % get projections
% [ p,mags ] = project( line1,v(1,:) );
% averages = centerAverages( averages );
% 
% figure;
% global hexagon_lattice_constant
% hold on;
% title('Disk 1')
% xlabel('Line 1 projection magnitude');
% % subplot(2,3,2);
% scatter(mags,averages(:,1,1)./max(averages(:,1,1)));
% myfun2 = @(y) cos(y*(pi/2)/(hexagon_lattice_constant/2)).^2;
% xlabel('Lattice plane offset');
% ylabel('Blinking intensity, disk 1');
% myfun = @(y) 0.85*abs(cos(y*(pi/2)/(hexagon_lattice_constant/2)));
% funvals = myfun(mags);
% plot(mags,funvals,'-r');
% myfun2 = @(y) 0.85*cos(y*(pi/2)/(hexagon_lattice_constant/2)).^2;
% funvals2 = myfun2(mags);
% plot(mags,funvals2,'-k');
% legend('Multislice Points','0.85*abs(cos)','cos^2');
% 
