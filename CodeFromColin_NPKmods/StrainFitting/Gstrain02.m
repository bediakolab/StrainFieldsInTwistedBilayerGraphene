function [sFit] = Gstrain02(sFit)

% Define graphene bilayer lattice


or = [258 233];
uv1 = [ ...
    96 32;
    16 104];
% tt = -1.2*pi/180;
% m = [cos(tt) -sin(tt);sin(tt) cos(tt)];
% uv1*m
uv2 = [ ...
    96.3 34;
    13.8 104.3];
basis = [
    1 1;
    -1 -1;
    1 -2;
    -1 2;
    -2 1;
    ];
% basis = [ ...
%     0 1;
%     0 -1;
%     -1 0;
%     -1 1;
%     1 -1;
%     ...
%     1 1;
%     -1 -1;
%     1 -2;
%     -1 2;
%     -2 1;
%     ];
   

intRange = [0 1000];


% Construct lattices
sFit.lat1 = [or; uv1];
sFit.lat2 = [or; uv2];
sFit.basis = [ones(size(basis,1),1) basis];
sFit.xyInit1 = sFit.basis*sFit.lat1;
sFit.xyInit2 = sFit.basis*sFit.lat2;

sFit.xyFit1 = sFit.xyInit1;
sFit.xyFit2 = sFit.xyInit2;


rPlot = 9;
t = linspace(0,2*pi,180+1);
ct = cos(t) * rPlot;
st = sin(t) * rPlot;


figure(11)
clf
imagesc(sFit.CBEDmax)
hold on
for a0 = 1:size(basis,1)
    plot(st + sFit.xyInit1(a0,2),...
        ct + sFit.xyInit1(a0,1),...
        'linewidth',1,'color',[1 0 0])
    plot(st + sFit.xyInit2(a0,2),...
        ct + sFit.xyInit2(a0,1),...
        'linewidth',1,'color',[0 0.7 1])
    
end
hold off
axis equal off
colormap(gray(256))
set(gca,'position',[0 0 1 1])
caxis(intRange)
% xlim([-50 50]+100)
% ylim([-50 50]+150)


end