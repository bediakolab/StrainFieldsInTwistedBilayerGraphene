% Hopefully this is the right way to think about it
function amppred = moirePredFun3(x0s,y0s,xor,yor,thetam)
thetam = deg2rad(thetam);
amppred = (x0s-xor).^2.*(2-2*cos(thetam)) + (y0s-yor).^2.*(2+2*cos(thetam)) + 2*(x0s-xor).*(y0s-yor).*(2*sin(thetam));
xcoords = (x0s-xor);
ycoords = (y0s-yor);
newcoords = [xcoords,ycoords];
rotmat = [cos(thetam),-sin(thetam);sin(thetam),cos(thetam)];
rotcoords = (rotmat*newcoords')';
uxpred = rotcoords(:,1) - newcoords(:,1);
uypred = rotcoords(:,2) - newcoords(:,2);
% uxpred = (x0s-xor)*(cos(thetam)-1) - (y0s-yor)*sin(thetam);
% uypred = (x0s-xor)*sin(thetam) + (y0s-yor)*(cos(thetam)+1);
amppred = uxpred.^2 + uypred.^2;
% shaperes = reshape(uxpred,[sqrt(numel(x0s)),sqrt(numel(x0s))]);
% figure;
% pcolor(shaperes);
% shading flat
end
