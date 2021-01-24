% Not sure why separable linear least squares isn't working
function residuals = moireFitFun2(ux,uy,x0s,y0s,xor,yor,thetam)
thetam = deg2rad(thetam);
uxpred = (x0s-xor)*(cos(thetam)-1) - (y0s-yor)*sin(thetam);
uypred = (x0s-xor)*sin(thetam) + (y0s-yor)*(cos(thetam)+1);
residuals = vertcat(ux,uy) - vertcat(uxpred,uypred);
shaperes = reshape(uxpred,[sqrt(numel(x0s)),sqrt(numel(x0s))]);
figure;
pcolor(shaperes);
shading flat
end
