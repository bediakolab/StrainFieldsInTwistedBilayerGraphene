function [posRefine] = Gfitting(sFit,stack4D,pos)

% Fitting of 2-disk overlaps
threshXray = 80;
flagPlot = false;
flagProgress = true;

% Fitting options
options = optimset( ...
    'TolX',1e-4,...
    'TolFun',1e-4,...
    'MaxIter',1e2,...
    'MaxFunEvals',1e3,...
    'display','off',...
    'UseParallel',0);   % NPK modified 02/03/2020 because NPK does not have the parallel computing toolbox.
%     'UseParallel',1); %,...
%     'ScaleProblem','jacobian');

% Coordinates
[ya,xa] = meshgrid(1:sFit.stackSize(2),1:sFit.stackSize(1));
CBEDtemplate = sFit.CBEDtemplateAmp( ...
    [(1:sFit.fitRadius) ((1-sFit.fitRadius):0)+sFit.stackSize(1)],...
    [(1:sFit.fitRadius) ((1-sFit.fitRadius):0)+sFit.stackSize(2)]);
CBEDtemplateFFT = fft2(CBEDtemplate);

% Main loop
Np = size(pos,1);
posRefine = zeros(Np,8,sFit.stackSize(3),sFit.stackSize(4));
for a0 = 1:Np
    % coordinates
    xy1 = pos(a0,1:2);
    xy2 = pos(a0,3:4);
    xv = round(mean([xy1(1) xy2(1)])) + sFit.fitVec;
    yv = round(mean([xy1(2) xy2(2)])) + sFit.fitVec;
    xCrop = xa(xv,yv);
    yCrop = ya(xv,yv);
    imageCrop = sFit.CBEDmax(xv,yv);
    imageCrop = filtXray(imageCrop,threshXray);
    
    % Make basis function
    %     basis = [xCrop(:) yCrop(:) CBEDtemplateFFT(:)];
    
    % Generate initial parameters
    intBG = median([imageCrop(:,1); imageCrop(:,end);
        imageCrop(:,1); imageCrop(:,end)]);
    intMax = max(imageCrop(:));
    % NPK modified 02/03/2020 to encase in double()
    a1 = sqrt(double(intMax - intBG)) * 0.5;
    a2 = a1 / sqrt(2);
    b2 = a1 / sqrt(2);
    coefsInit = [ ...
        xy1 xy2 intBG a1 a2 b2];
    
    % Bounds
    lb = [coefsInit(1:4) - sFit.fitDeltaXY ...
        0 -Inf -Inf -Inf];
    ub = [coefsInit(1:4) + sFit.fitDeltaXY ...
        Inf Inf Inf Inf];
    
    % Testing
    %     imageOut = drawCBED(coefsInit,basis);
    
    % Loop over probe positions
    for ax = 1:sFit.stackSize(3)
        fprintf('Optimizing data row %d of %d.\n',ax,sFit.stackSize(3));
        for ay = 1:sFit.stackSize(4)
            imageCrop = stack4D(xv,yv,ax,ay);
            imageCrop = filtXray(imageCrop,threshXray);
            basis = [xCrop(:) yCrop(:) CBEDtemplateFFT(:)];
            coefs = coefsInit;
            
            % Fitting
            try
                % NPK modified 02/03/2020 to encase everything in double(),
                % avoiding errors from lsqcurvefit
                coefs = lsqcurvefit(@drawCBED,double(coefs),basis,double(imageCrop),double(lb),double(ub),options);
            catch err
                disp(['Peak ' num2str(a0) ...
                    ', ax = ' num2str(ax) ...
                    ', ay = ' num2str(ay) ...
                    ', err ' err])
                coefs(5:end) = 0;
            end
            
            % Output
            posRefine(a0,:,ax,ay) = coefs;
            
            % plotting
            if flagPlot == true
                imageOut = drawCBED(coefs,basis);
                figure(45)
                clf
                imagesc([imageCrop imageOut])
                axis equal off
                colormap(jet(256))
                caxis([min(imageCrop(:)) max(imageCrop(:))])
                drawnow
            end
            
            if flagProgress == true
                comp = ((ay / sFit.stackSize(4) ...
                    + ax - 1) / sFit.stackSize(3) ...
                    + a0 - 1) / Np;
                progressbar(comp,2);
            end
        end
    end
    
    
end

end


function [imageOut] = drawCBED(coefs,basis)
% Probe generation function

N = [1 1]*sqrt(size(basis,1));
imageOut = zeros(N);

x = coefs(1) - basis(1,1);
y = coefs(2) - basis(1,2);
xF = floor(x);
yF = floor(y);
dx = x - xF;
dy = y - yF;
inds = sub2ind(N,[xF xF+1 xF xF+1],[yF yF yF+1 yF+1]);
imageOut(inds) = imageOut(inds) + ...
    coefs(6) * [(1-dx)*(1-dy) dx*(1-dy) (1-dx)*dy dx*dy];

x = coefs(3) - basis(1,1);
y = coefs(4) - basis(1,2);
xF = floor(x);
yF = floor(y);
dx = x - xF;
dy = y - yF;
inds = sub2ind(N,[xF xF+1 xF xF+1],[yF yF yF+1 yF+1]);
imageOut(inds) = imageOut(inds) + ...
    (coefs(7) + 1i*coefs(8)) * [(1-dx)*(1-dy) dx*(1-dy) (1-dx)*dy dx*dy];

imageOut(:) = coefs(5) + ...
    abs(ifft2(fft2(imageOut) .* reshape(basis(:,3),N))).^2;

end