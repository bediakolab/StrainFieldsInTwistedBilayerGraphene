function [sFit] = Gstrain06(sFit)


strainRange = [-1 1]*0.008;

% Generate strain maps for bilayer graphene

distWeight = 5 - 2;

% Ref lattices
N = sFit.stackSize;
basis = sFit.basis;
Np = size(basis,1);
xy1 = sFit.posRefine(:,1:2,:,:);
xy2 = sFit.posRefine(:,3:4,:,:);


x1 = reshape(xy1(:,1,:,:),[Np N(3)*N(4)]);
y1 = reshape(xy1(:,2,:,:),[Np N(3)*N(4)]);
x2 = reshape(xy2(:,1,:,:),[Np N(3)*N(4)]);
y2 = reshape(xy2(:,2,:,:),[Np N(3)*N(4)]);
x1ref = median(x1,2);
y1ref = median(y1,2);
x2ref = median(x2,2);
y2ref = median(y2,2);

lat1ref = basis \ [x1ref y1ref];
lat2ref = basis \ [x2ref y2ref];
or = (lat1ref(1,:) + lat2ref(1,:))/2;
uv1ref = basis(:,2:3) \ [x1ref-or(1) y1ref-or(2)];
uv2ref = basis(:,2:3) \ [x2ref-or(1) y2ref-or(2)];
sFit.strainLat1 = [or; uv1ref];
sFit.strainLat2 = [or; uv2ref];


% compute strain maps
sFit.strainMaps1 = zeros(N(3),N(4),8);
sFit.strainMaps2 = zeros(N(3),N(4),8);
for ax = 1:N(3)
    for ay = 1:N(4)     
        xy1fit = xy1(:,1:2,ax,ay);
        xy2fit = xy2(:,1:2,ax,ay);
        %         w1 = ones(Np,1);
        %         w2 = ones(Np,1);
        dist1_2 = ( ...
            (xy1fit(:,1) - x1ref).^2 + ...
            (xy1fit(:,2) - y1ref).^2);
        dist2_2 = ( ...
            (xy2fit(:,1) - x2ref).^2 + ...
            (xy2fit(:,2) - y2ref).^2);
        w1 = exp(dist1_2/(-2*distWeight^2));
        w2 = exp(dist2_2/(-2*distWeight^2));
        
        % Fit lattices
        uv1 = lscov(basis(:,2:3),xy1fit-or,w1);
        uv2 = lscov(basis(:,2:3),xy2fit-or,w2);
        
        % Measure strains in lat 1
        m1 = uv1ref \ uv1;        
        exx = 1 - m1(1,1);
        eyy = 1 - m1(2,2);
        exy = (-0.5)*(m1(1,2) + m1(2,1));
        theta = 0.5*(m1(2,1) - m1(1,2));
        sFit.strainMaps1(ax,ay,1:4) = [exx eyy exy theta];
        sFit.strainMaps1(ax,ay,5:8) = m1(:);
        
        % Measure strains in lat 2
        m2 = uv2ref \ uv2;
        exx = 1 - m2(1,1);
        eyy = 1 - m2(2,2);
        exy = (-0.5)*(m2(1,2) + m2(2,1));
        theta = 0.5*(m2(2,1) - m2(1,2));
        sFit.strainMaps2(ax,ay,1:4) = [exx eyy exy theta];
        sFit.strainMaps2(ax,ay,5:8) = m2(:);
    end
    
    comp = ax / N(3);
    progressbar(comp,2);
end

b = 5;
mf = [1 1]+3;
Ip1 = [ ...
    medfilt2(sFit.strainMaps1(:,:,1),mf,'symmetric') ...
    ones(N(3),b) ...
    medfilt2(sFit.strainMaps1(:,:,2),mf,'symmetric') ...
    ones(N(3),b) ...
    medfilt2(sFit.strainMaps1(:,:,3),mf,'symmetric') ...
    ones(N(3),b) ...
    medfilt2(sFit.strainMaps1(:,:,4),mf,'symmetric')];
Ip2 = [ ...
    medfilt2(sFit.strainMaps2(:,:,1),mf,'symmetric') ...
    ones(N(3),b) ...
    medfilt2(sFit.strainMaps2(:,:,2),mf,'symmetric') ...
    ones(N(3),b) ...
    medfilt2(sFit.strainMaps2(:,:,3),mf,'symmetric') ...
    ones(N(3),b) ...
    medfilt2(sFit.strainMaps2(:,:,4),mf,'symmetric')];
    

figure(1)
clf
imagesc([Ip1 - Ip2]*100)
axis equal off
colormap(jet(356))
caxis(strainRange*100)
colorbar


% u1 = lat1ref(3,:)
% u2 = lat2ref(3,:)
% acos(sum(u1.*u2)/norm(u1)/norm(u2))*180/pi;

% Plotting
figure(11)
clf
hold on
x = xy1(:,1,:,:);
y = xy1(:,2,:,:);
scatter(y(:),x(:),'k.')

x = xy2(:,1,:,:);
y = xy2(:,2,:,:);
scatter(y(:),x(:),'b.')

hold off

axis equal 
set(gca,'ydir','reverse')

end