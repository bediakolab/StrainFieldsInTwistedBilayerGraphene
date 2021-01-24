function [] = Gstrain07(sFit)


% Plotting function for strains

flagLines = 1;

p1 = [11 33];
p2 = [30 33];
p3 = [21 22];
p4 = [2 22];
p5 = [39 22];
p6 = [11 12];
p7 = [30 12];

L = [ ...
    p1 p2;
    p2 p3;
    p3 p1; 
    ...
    p1 p4;
    p4 p3;
    ...
    p5 p2;
    p5 p3;
    ...
    p6 p4;
    p6 p3;
    ...
    p7 p3;
    p7 p5;
    p7 p6];

strainRange = [-1 1]*0.01;
thetaRange = [-1 1]*0.5 * (pi/180);
filtSize = [1 1]*3;

% Get data
exx1 = sFit.strainMaps1(:,:,1);
eyy1 = sFit.strainMaps1(:,:,2);
exy1 = sFit.strainMaps1(:,:,3);
the1 = sFit.strainMaps1(:,:,4);
exx2 = sFit.strainMaps2(:,:,1);
eyy2 = sFit.strainMaps2(:,:,2);
exy2 = sFit.strainMaps2(:,:,3);
the2 = sFit.strainMaps2(:,:,4);

% Filtering
exx1 = medfilt2(exx1,filtSize,'symmetric');
eyy1 = medfilt2(eyy1,filtSize,'symmetric');
exy1 = medfilt2(exy1,filtSize,'symmetric');
the1 = medfilt2(the1,filtSize,'symmetric');
exx2 = medfilt2(exx2,filtSize,'symmetric');
eyy2 = medfilt2(eyy2,filtSize,'symmetric');
exy2 = medfilt2(exy2,filtSize,'symmetric');
the2 = medfilt2(the2,filtSize,'symmetric');





h = figure(11);
clf
set(gcf,'color','w','outerposition',[510 120 576+512 512+384])

b = [0.02 0.02];

ha = cell(3,4);

ha{1,1} = axes('position',[0/4+b(1) 0.6+b(2) 1/4-2*b(1) 0.3-2*b(2)]);
ha{1,2} = axes('position',[1/4+b(1) 0.6+b(2) 1/4-2*b(1) 0.3-2*b(2)]);
ha{1,3} = axes('position',[2/4+b(1) 0.6+b(2) 1/4-2*b(1) 0.3-2*b(2)]);
ha{1,4} = axes('position',[3/4+b(1) 0.6+b(2) 1/4-2*b(1) 0.3-2*b(2)]);

ha{2,1} = axes('position',[0/4+b(1) 0.3+b(2) 1/4-2*b(1) 0.3-2*b(2)]);
ha{2,2} = axes('position',[1/4+b(1) 0.3+b(2) 1/4-2*b(1) 0.3-2*b(2)]);
ha{2,3} = axes('position',[2/4+b(1) 0.3+b(2) 1/4-2*b(1) 0.3-2*b(2)]);
ha{2,4} = axes('position',[3/4+b(1) 0.3+b(2) 1/4-2*b(1) 0.3-2*b(2)]);

ha{3,1} = axes('position',[0/4+b(1) 0.0+b(2) 1/4-2*b(1) 0.3-2*b(2)]);
ha{3,2} = axes('position',[1/4+b(1) 0.0+b(2) 1/4-2*b(1) 0.3-2*b(2)]);
ha{3,3} = axes('position',[2/4+b(1) 0.0+b(2) 1/4-2*b(1) 0.3-2*b(2)]);
ha{3,4} = axes('position',[3/4+b(1) 0.0+b(2) 1/4-2*b(1) 0.3-2*b(2)]);


for ax = 1:3
    switch ax
        case 1
            I1 = exx1;
            I2 = eyy1;
            I3 = exy1;
            I4 = the2;
        case 2
            I1 = exx2;
            I2 = eyy2;
            I3 = exy2;
            I4 = the2;
        case 3
            I1 = exx1 - exx2;
            I2 = eyy1 - eyy2;
            I3 = exy1 - exy2;
            I4 = the1 - the2;
    end
    
    
    for ay = 1:4
        set(h,'currentaxes',ha{ax,ay});
        
        switch ay
            case 1
                I = I1;
                cr = strainRange;
            case 2
                I = I2;
                cr = strainRange;
            case 3
                I = I3;
                cr = strainRange;
            case 4
                I = I4;     
                cr = thetaRange;
        end
        
        imagesc(I - mean(I(:)))
        hold on
        if flagLines == true
            for a0 = 1:size(L,1)
                line(L(a0,[2 4]),L(a0,[1 3]),...
                    'linewidth',2,'color','k')
            end
        end
        hold off
        axis equal off
        caxis(cr);
        

        
        
    end
end




colormap(jet(245))





end