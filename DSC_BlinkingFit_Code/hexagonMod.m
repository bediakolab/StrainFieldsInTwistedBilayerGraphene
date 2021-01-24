function [ output_args ] = hexagonMod( v1,v2,data )
% Nathanael Kazmierczak
% Built for the case (as with DSC lattice) where v1 is straight up, v2 is
% at an angle to the right.

[r,c] = size(v1);
if r > c
    v1 = v1';
    v2 = v2';
end

h = zeros(6,2);
h(1,:) = v1;
h(2,:) = v2;
h(3,:) = v2-v1;
h(4,:) = -v1;
h(5,:) = -v2;
h(6,:) = v1-v2;

vl1 = vectorLine(v2-v1,v1);
data(vl1 < data) = data(vl1 < data) - (v1+v2);



end

