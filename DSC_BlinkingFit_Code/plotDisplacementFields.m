function plotDisplacementFields( processed_dfield, see_options )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

dampf = (processed_dfield(:,:,1).^2 + processed_dfield(:,:,2).^2).^(1/2);
danglef = atan2(processed_dfield(:,:,2),processed_dfield(:,:,1));
figh = figure;
subplot(2,3,1)
imagesc(dampf); colormap(fire); axis equal; colorbar;
title('Reduced processed amp');
set(gca,'ydir','normal');
subplot(2,3,2)
imagesc(danglef); colormap(hsv); axis equal; colorbar;
title('Reduced processed angle');
set(gca,'ydir','normal');

subplot(2,3,4)
pd1 = processed_dfield(:,:,1);
pd2 = processed_dfield(:,:,2);
scatter( pd1(:), pd2(:), 'filled','r' );
axh = gca;
plotFullDisplacementHexagons( axh );

subplot(2,3,5)
xbase1 = 1:size(processed_dfield,1);
ybase1 = 1:size(processed_dfield,2);
[xspace1,yspace1] = meshgrid(xbase1,ybase1);
quiver(xspace1(:),yspace1(:),pd1(:),pd2(:),0);
xlim([1,xbase1(end)]);
ylim([1,ybase1(end)]);

input('Please zoom to a vector and press any key when done.');
figure(figh);
disp('Please click once on a vector');
[xinit,yinit] = ginput(1);
x = round(yinit);
y = round(xinit);
pd1val = pd1(x,y);
pd2val = pd2(x,y);
hold on;
quiver(xspace1(x,y),yspace1(x,y),pd1(x,y),pd2(x,y),0,'r');
figure(figh);







% set(figh,'Position',[0,0,1000,1000])
if see_options
    position = [pd1val,pd2val];
    
    int_bound = 6;
    [ v1, v2, hexagon_lattice_constant ] = getDSCBasisVectors();
%     axes(plothandle);
    bound_handling_flag = 1;
    scaling_constant = 1;
    [ pred_vals ] = trigFittingFunctions( position, scaling_constant, bound_handling_flag );
    
    positions = [v1;v2;v2-v1;-v1;-v2;v1-v2;...
        2*v1-v2;v1+v2;2*v2-v1;v2-2*v1;-v1-v2;v1-2*v2];
    n = 500;
    xbase = linspace(-int_bound,int_bound,n);
    ybase = linspace(-int_bound,int_bound,n);
    [xspace,yspace] = meshgrid(xbase,ybase);
    r = 0.2;
    
    basis = permn(-int_bound:int_bound,2);
    basis(basis(:,1) == 0 & basis(:,2) == 0,:) = [];
    % basis = [2 2
    %          2 1
    %          2 0
    %          2 -1
    %          2 -2
    %          1 2
    %          1 1
    %          1 0
    %          1 -1
    %          1 -2
    %          0 2
    %          0 1
    %
    %          0 -1
    %          0 -2
    %          -1 2
    %          -1 1
    %          -1 0
    %          -1 -1
    %          -1 -2
    %          -2 2
    %          -2 1
    %          -2 0
    %          -2 -1
    %          -2 -2];
    
    w1 = v1 + v2;
    w2 = 2*v2 - v1;
    W = [w1',w2'];
    genpoints = W*basis';
    genpoints_plus = genpoints + repmat(position',1,size(genpoints,2));
    genpoints_minus = genpoints - repmat(position',1,size(genpoints,2));
    genpoints_all = [genpoints_plus,genpoints_minus,-position']';
    
    subplot(2,3,6);
    plot(genpoints_all(:,1),genpoints_all(:,2),'ro');
    hold on
    plot(position(1),position(2),'bo');
    % ch = viscircles(position,RADIUS);
    xlim([-4,4]);
    ylim([-4,4]);
    axh = gca;
    plotFullDisplacementHexagons(axh);
    
    while true
        [a,b] = ginput(1);
        [~,idx] = min(sqrt((genpoints_all(:,1) - a).^2 + (genpoints_all(:,2) - b).^2));
        newvec = genpoints_all(idx,:);
        
        
        % plot on new axes
        pd1_copy = pd1;
        pd2_copy = pd2;
        pd1_copy(x,y) = newvec(1);
        pd2_copy(x,y) = newvec(2);
        subplot(2,3,3);
        quiver(xspace1(:),yspace1(:),pd1_copy(:),pd2_copy(:),0,'b');
        hold on
        quiver(xspace1(x,y),yspace1(x,y),pd1_copy(x,y),pd2_copy(x,y),0,'r');
        tf = input('Would you like to try another equivalent vector? 1/0 ');
        if ~tf
            break
        end
            
    end
end


end

