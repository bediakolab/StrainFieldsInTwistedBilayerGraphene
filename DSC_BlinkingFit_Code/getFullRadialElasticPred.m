function [predmat] = getFullRadialElasticPred(c,matsize,function_string,const_vals)
% Inner function for computing model values for elastic scattering for
% multiple fitting functions.
%
% Nathanael Kazmierczak, 03/12/2020
if nargin < 3
    const_vals = [];
end

switch function_string
    case 'power law'
        % c(1) is x0, c(2) is y0, c(3) is log10(A), and c(4) is B (the exponent).
        A = 10^c(3);
        B = c(4);
        x0 = c(1);
        y0 = c(2);
        xbase = 1:matsize(2);
        ybase = 1:matsize(1);
        [xspace,yspace] = meshgrid(xbase,ybase);
        xspace_centered = xspace - x0;
        yspace_centered = yspace - y0;
        radii = sqrt(xspace_centered.^2 + yspace_centered.^2);
        predmat = A*radii.^B;
    case 'lorentzian'
        % Assumes a spacing of one in each pixel direction.
        % c(1) is x0, c(2) is y0, c(3) is A (scaling), c(4) is B (the width parameter).
        A = c(3);
        B = c(4);
        x0 = c(1);
        y0 = c(2);
        xbase = 1:matsize(2);
        ybase = 1:matsize(1);
        [xspace,yspace] = meshgrid(xbase,ybase);
        xspace_centered = xspace - x0;
        yspace_centered = yspace - y0;
        radii = sqrt(xspace_centered.^2 + yspace_centered.^2);
        predmat = A*(1/pi)*(0.5*B)./(radii.^2 + (0.5*B)^2);
    case 'gaussian'
        A = c(3);
        fG = c(4);
        x0 = c(1);
        y0 = c(2);
        xbase = 1:matsize(2);
        ybase = 1:matsize(1);
        [xspace,yspace] = meshgrid(xbase,ybase);
        xspace_centered = xspace - x0;
        yspace_centered = yspace - y0;
        radii = sqrt(xspace_centered.^2 + yspace_centered.^2);
        sigma = fG/(2*sqrt(2*log(2)));  % Gaussian width parameter
        predmat = A*(1/(sigma*sqrt(2*pi)))*exp(-0.5*(radii./sigma).^2);
    case 'pseudo-voigt'
        % Assumes a spacing of one in each pixel direction. c(1) is x0, c(2) is y0,
        % c(3) is A (scaling), c(4) is fG, c(5) is fL (the full width at half max
        % for Gaussian and Lorentzian, respectively).
        A = c(3);
        fG = c(4);
        fL = c(5);
        x0 = c(1);
        y0 = c(2);
        xbase = 1:matsize(2);
        ybase = 1:matsize(1);
        [xspace,yspace] = meshgrid(xbase,ybase);
        xspace_centered = xspace - x0;
        yspace_centered = yspace - y0;
        radii = sqrt(xspace_centered.^2 + yspace_centered.^2);
        sigma = fG/(2*sqrt(2*log(2)));  % Gaussian width parameter
        gamma = fL/2;  % Lorentz width parameter
        f = (fG^5 + 2.69269*fG^4*fL + 2.42843*fG^3*fL^2 + 4.47163*fG^2*fL^3 + 0.07842*fG*fL^4 + fL^5)^1/5;
        eta = 1.36603*(fL/f) - 0.47719*(fL/f)^2 + 0.11116*(fL/f)^3;
        lorentz_predmat = A*(1/pi)*(0.5*gamma)./(radii.^2 + (0.5*gamma)^2);
        gaussian_predmat = A*(1/(sigma*sqrt(2*pi)))*exp(-0.5*(radii./sigma).^2);
        predmat = eta*lorentz_predmat + (1-eta)*gaussian_predmat;
    case 'pseudo-voigt unconstrained'
        % Assumes a spacing of one in each pixel direction. c(1) is x0, c(2) is y0,
        % c(3) is A (scaling), c(4) is fG, c(5) is fL (the full width at half max
        % for Gaussian and Lorentzian, respectively).
        % This function is unconstrained because c(6) is the linear
        % combination weighting parameter eta, which otherwise can be computed
        % from the PV expression.
        AG = c(3);
        AL = c(4);
        fG = c(5);
        fL = c(6);
        x0 = c(1);
        y0 = c(2);
        xbase = 1:matsize(2);
        ybase = 1:matsize(1);
        [xspace,yspace] = meshgrid(xbase,ybase);
        xspace_centered = xspace - x0;
        yspace_centered = yspace - y0;
        radii = sqrt(xspace_centered.^2 + yspace_centered.^2);
        sigma = fG/(2*sqrt(2*log(2)));  % Gaussian width parameter
        gamma = fL/2;  % Lorentz width parameter
        % f = (fG^5 + 2.69269*fG^4*fL + 2.42843*fG^3*fL^2 + 4.47163*fG^2*fL^3 + 0.07842*fG*fL^4 + fL^5)^1/5;
        % eta = 1.36603*(fL/f) - 0.47719*(fL/f)^2 + 0.11116*(fL/f)^3;
        lorentz_predmat = AL*(1/pi)*(0.5*gamma)./(radii.^2 + (0.5*gamma)^2);
        gaussian_predmat = AG*(1/(sigma*sqrt(2*pi)))*exp(-0.5*(radii./sigma).^2);
        predmat = lorentz_predmat + gaussian_predmat;
    case 'lorentz + sinc'
    case 'ringed gaussian'
        I0 = c(1);
        Ipeak = c(2);
        Iring = c(3);
        rRing = c(4);
        fwhm = c(5);
        fwhm1 = c(6);
        fwhm2 = c(7);
        x0 = const_vals(1);
        y0 = const_vals(2);
        xbase = 1:matsize(2);
        ybase = 1:matsize(1);
        [xspace,yspace] = meshgrid(xbase,ybase);
        xspace_centered = xspace - x0;
        yspace_centered = yspace - y0;
        r = sqrt(xspace_centered.^2 + yspace_centered.^2);
        s = fwhm/(2*sqrt(2*log(2)));  % Gaussian width parameter
        s1 = fwhm1/(2*sqrt(2*log(2)));
        s2 = fwhm2/(2*sqrt(2*log(2)));
        
        predmat = I0 + Ipeak.*exp(-0.5*(r./s).^2) + ...
            (r<rRing).*Iring.*exp(-0.5*((r-rRing)./s1).^2) + ...
            (r>=rRing).*Iring.*exp(-0.5*((r-rRing)./s2).^2);
            
end



end

