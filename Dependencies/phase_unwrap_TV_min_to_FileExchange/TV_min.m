%% *Function to minimize noise on input image or 2-dimensional data matrix by using TV (Total Variation) minimization algorithm.[1]*
% All equations here are from Ref.[2], although Ref.[2] is application paper of
% TV (Total Variation) minimization algorithm originally suggested by
% Ref.[1].
% How to use: denoised_image = TV_minimization(varargin) where varargin =
% input_img or input_image,mu. input_img is 2-dimensional image data matrix
% (input_img = imread(your image) and mu is regularization parameter
% influencing noise-reduction process. Ref.[2] claims that mu between 0.5
% to 5 is good. Without mu as a argument of the function, mu is set to be 2
% as a default value.
function denoised_image = TV_min(varargin)
tau = 0.25; 
input_img = varargin{1}; 
px_1th = zeros(size(input_img)); py_1th = zeros(size(input_img)); % x and y-component of zero vector field as initiation parameter of iternation (Eq.(9)).
if (length(varargin) == 2)
    mu = varargin{2};
else
    mu = 2;
end
% 1st-iteration given by Eq.(9)
[com_factor_x, com_factor_y] = gradient(divergence(px_1th,py_1th) - input_img/mu); % common factor to both numerator and denominator of Eq.(9).
norm_com_factor = (com_factor_x.^2 + com_factor_y.^2).^(1/2); % norm of the common factor in denominator.
px_nth = (px_1th + tau*com_factor_x)./(1 + tau.*norm_com_factor); py_nth = (py_1th + tau*com_factor_y)./(1 + tau.*norm_com_factor);
% Rest iteration
while (1)
    [com_factor_x, com_factor_y] = gradient(divergence(px_nth,py_nth) - input_img/mu); % common factor to both numerator and denominator of Eq.(9).
    norm_com_factor = (com_factor_x.^2 + com_factor_y.^2).^(1/2); % norm of the common factor in denominator.
    px_n_plus_1th = (px_nth + tau*com_factor_x)./(1 + tau.*norm_com_factor); py_n_plus_1th = (py_nth + tau*com_factor_y)./(1 + tau.*norm_com_factor);
    if (max(max((divergence(px_n_plus_1th,py_n_plus_1th) - divergence(px_nth,py_nth)))) < 0.002) % Iteration-stopping criteria mentioned in p.7 of Ref.[1].
        break;
    else
        px_nth = px_n_plus_1th; py_nth = py_n_plus_1th;
    end
end
pi_muk = mu*(divergence(px_n_plus_1th,py_n_plus_1th)); % nonlinear projectin calculated by using Eq.
denoised_image = input_img - pi_muk;
end
%% Reference
% # ANTONIN CHAMBOLLE, ¡°An Algorithm for Total Variation Minimization and Applications,¡± J. Math. Imaging Vis., vol. 20, no. 1/2, pp. 89?97, Jan. 2004.
% # H. Y. H. Huang, L. Tian, Z. Zhang, Y. Liu, Z. Chen, and G. Barbastathis, ¡°Path-independent phase unwrapping using phase gradient and total-variation (TV) denoising,¡± Opt. Express, vol. 20, no. 13, p. 14075, Jun. 2012.