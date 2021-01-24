%% *Function to convert unwrapped phase into wrapped phase (Eq.(2) of Ref.[1])*
%% Function definition
function wrapped_phase = phase_wrap(unwrapped_phase)
[m, n] = size(unwrapped_phase);
wrapped_phase = 2*atan(sin(unwrapped_phase)./(1 + cos(unwrapped_phase)));
for i = 1 : m
    for j = 1 : n
        if (cos(unwrapped_phase(i,j)) == -1)
            wrapped_phase(i,j) = pi;
        end
    end
end
%% Reference
% # H. Y. H. Huang, L. Tian, Z. Zhang, Y. Liu, Z. Chen, and G. Barbastathis, ¡°Path-independent phase unwrapping using phase gradient and total-variation (TV) denoising,¡± Opt. Express, vol. 20, no. 13, p. 14075, Jun. 2012.