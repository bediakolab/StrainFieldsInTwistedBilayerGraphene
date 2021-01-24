function makeProjectedBlinkingPlot( disk_num, averages, DSC_values)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

if numel(disk_num > 1)
    switch numel(disk_num)
        case 6
            r = 2;
            c = 3;
    end
end

[ ~, ~, hexagon_lattice_constant ] = getDSCBasisVectors();
[ dict_function, v ] = getLatticePlaneNormalVectors();


figure
for i = 1:numel(disk_num)
    this_disk_num = disk_num(i);
    
    if numel(disk_num > 1)
        subplot(r,c,i);
    end
    
    
    if this_disk_num < 7
        % Inner
        %     plane_distance = hexagon_lattice_constant*3/2;
        plane_distance = 3*hexagon_lattice_constant/2/sqrt(3);
    else
        % Outer
        %     plane_distance = hexagon_lattice_constant*sqrt(3)/2;
        plane_distance = hexagon_lattice_constant/2;
    end
    
    % get projections
    [ p,mags ] = project( DSC_values,v(this_disk_num,:) );
    averages = centerAverages( averages );
    
    mul = (pi/2)./(plane_distance/2);
    
    
    hold on;
    title(sprintf('Disk %d',this_disk_num))
    xlabel('Projection magnitude');
    % subplot(2,3,2);
    scatter(mags,averages(:,this_disk_num)./max(averages(:,this_disk_num)));
    % myfun2 = @(y) cos(y*mul).^2;
    xlabel('Lattice plane offset');
    ylabel('Blinking intensity');
    % myfun = @(y) 0.9*abs(cos(y*mul));
    % funvals = myfun(mags);
    % plot(mags,funvals,'-r');
    scaling = 0.9;
    myfun2 = @(y) scaling*cos(y*mul).^2;
    funvals2 = myfun2(mags);
    plot(mags,funvals2,'-k');
    legend('Multislice Points',sprintf('%.2f*cos^2(offset)',scaling));
    
    
end

end

