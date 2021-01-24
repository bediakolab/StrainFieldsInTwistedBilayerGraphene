function [sFit] = Gstrain04NPKmod(sFit,stack4D)


sFit.fitRadius = 16;
sFit.fitVec = (1-sFit.fitRadius):sFit.fitRadius;
sFit.fitDeltaXY = 5;

if nargout == 0
    indPlot = 1 + 5;
    x = round(sFit.xyFit1(indPlot,1)) + sFit.fitVec;
    y = round(sFit.xyFit1(indPlot,2)) + sFit.fitVec;
    % Ip = sFit.CBEDmax(x,y);
    Ip = stack4D(x,y,22,19);
    Ip = filtXray(Ip,80);
    figure(11)
    clf
    imagesc(Ip)
    axis equal off
    colormap(gray(256))
else
    %     pos = [sFit.xyInit1(7,:) sFit.xyInit2(7,:)];
    pos = [sFit.xyInit1 sFit.xyInit2];
    if isempty(sFit.posRefine)
        sFit.posRefine = Gfitting(sFit,stack4D,pos);
    else
        % Cat logic
        Gfit_result = Gfitting(sFit,stack4D,pos);
        sFit.posRefine = cat(3,sFit.posRefine,Gfit_result);
    end
end


end
