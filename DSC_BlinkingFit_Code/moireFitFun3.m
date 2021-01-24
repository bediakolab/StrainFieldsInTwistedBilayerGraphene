% Hopefully this is the right way to think about it
function residuals = moireFitFun3(usqamp,x0s,y0s,xor,yor,thetam)
% thetam = deg2rad(thetam);
% amppred = (x0s-xor).^2.*(2-2*cos(thetam)) + (y0s-yor).^2.*(2+2*cos(thetam)) + 2*(x0s-xor).*(y0s-yor).*(2*sin(thetam));
amppred = moirePredFun3(x0s,y0s,xor,yor,thetam);
residuals = usqamp - amppred;
% shaperes = reshape(uxpred,[sqrt(numel(x0s)),sqrt(numel(x0s))]);
% figure;
% pcolor(shaperes);
% shading flat
end
