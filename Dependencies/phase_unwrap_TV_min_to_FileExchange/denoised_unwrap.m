%% *Function to obtain denoised-unwrapped phase*
% Called functions in execution: phase_wrap, TV_min.
function denoised_unwrapped_phase = denoised_unwrap(wrapped_phase)
%% Calculating denoised derivatives (components of gradient) of unwrapped phase by using TV (Total Variation) minization algorithm [1][2]
gradient_x_wrapped_phase = diff(wrapped_phase,1,2); gradient_y_wrapped_phase = diff(wrapped_phase,1,1); 
gradient_x_unwrapped_phase = phase_wrap(gradient_x_wrapped_phase); gradient_y_unwrapped_phase = phase_wrap(gradient_y_wrapped_phase); % Implementation of Eq.(3) and (4) of Ref.[1].
%figure(11); subplot(2,1,1); imagesc(gradient_x_wrapped_phase); colormap(gray); colorbar;
subplot(2,1,2); imagesc(gradient_x_unwrapped_phase); colormap(gray); colorbar;
denoised_gradient_x_unwrapped_phase = TV_min(gradient_x_unwrapped_phase); denoised_gradient_y_unwrapped_phase = TV_min(gradient_y_unwrapped_phase); % Applying TV-minimization algorithm to derivatives.
%figure(12); subplot(2,1,1); surf(gradient_x_unwrapped_phase(:,:)); subplot(2,1,2); surf(denoised_gradient_x_unwrapped_phase(:,:));
%% Integration of denoised gradient to obtain denoised-unwrapped phase
gradient_x = denoised_gradient_x_unwrapped_phase; gradient_y = denoised_gradient_y_unwrapped_phase;
denoised_unwrapped_phase = zeros(size(wrapped_phase));
denoised_unwrapped_phase(1,1) = wrapped_phase(1,1);
for j = 2 : size(denoised_unwrapped_phase,2)
    denoised_unwrapped_phase(1,j) = denoised_unwrapped_phase(1,j-1) + gradient_x(1,j-1); % Integration along 1st row.
end
for i = 2 : size(denoised_unwrapped_phase,1)
    denoised_unwrapped_phase(i,1) = denoised_unwrapped_phase(i-1,1) + gradient_y(i-1,1); % Integration along 1st column.
end
for i = 2 : size(denoised_unwrapped_phase,1) % Integration for rest matrix body.
    for j = 2: size(denoised_unwrapped_phase,2)
        denoised_unwrapped_phase(i,j) = denoised_unwrapped_phase(i,j-1) + gradient_x(i,j-1);
        denoised_unwrapped_phase(i,j) = denoised_unwrapped_phase(i,j) + (denoised_unwrapped_phase(i-1,j) + gradient_y(i-1,j));
        denoised_unwrapped_phase(i,j) = denoised_unwrapped_phase(i,j)/2; % This mean-integration gives more reliable result than using derivative of only one direction.
    end
end
end
%% Reference
% # H. Y. H. Huang, L. Tian, Z. Zhang, Y. Liu, Z. Chen, and G. Barbastathis, ¡°Path-independent phase unwrapping using phase gradient and total-variation (TV) denoising,¡± Opt. Express, vol. 20, no. 13, p. 14075, Jun. 2012.
% # ANTONIN CHAMBOLLE, ¡°An Algorithm for Total Variation Minimization and Applications,¡± J. Math. Imaging Vis., vol. 20, no. 1/2, pp. 89?97, Jan. 2004.
