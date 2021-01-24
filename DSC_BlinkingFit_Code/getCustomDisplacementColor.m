function [ RGB_color_stack, HSV_color_stack, H_weights ] = getCustomDisplacementColor( displacement_field, SP_hsv, AB_hsv, SP_uniqueness_flag, gamma )
% Function for making 2D colormapping of a displacement field.
% H_weights is a variable for extracting the modular angular weights, to
% permit angle-based filtering in a rigorous (correct) way.
%
% Nathanael Kazmierczak, last modified 05/19/2020

% Original guesses 
% SP_hsv = [0.70833,0.5,1];
% AB_hsv = [0.56667,0.28,1];
% New from fiddling
if nargin < 2
    SP_hsv = [0.8,0.8,1];
    AB_hsv = [0.5,0.3,1];
end
if isempty(SP_uniqueness_flag)
    SP_uniqueness_flag = 0;
end
if nargin < 5
    gamma = 1.3;
end
if SP_uniqueness_flag
    %     AB_hsv = [0,  0,    1;
    %               0.33,  0,    1
    %               0.66,  0,    1];
    %     SP_hsv = [0,    1,      1;
    %               0.33, 1,      1;
    %               0.66, 1,      1];
    % NPK's original warmer hues color scheme.
    % % % h = [0.75,1,1.07];
    % % % AB_hsv = [h(1),0,1;
    % % %     h(2),0,1;
    % % %     h(3),0,1];
    % % % SP_hsv = [
    % % %     h(1),1,1;
    % % %     h(2),0.9,1;
    % % %     h(3),1,1];
    if isempty(SP_hsv)
        h = [0.75,1,1.07];
        AB_hsv = [h(1),0,1;
            h(2),0,1;
            h(3),0,1];
        SP_hsv = [
            h(1),1,1;
            h(2),0.9,1;
            h(3),1,1];
    end
end
% Testing
% SP_hsv = [0.8,0.8,1];
% AB_hsv = [0.5,0.3,1];

[ v1, v2, hlc ] = getDSCBasisVectors();
side_dimension = hlc/2;
vertex_dimension = hlc/sqrt(3);

tangents = displacement_field(:,:,2)./displacement_field(:,:,1);
tangents(isnan(tangents)) = inf;  % to get the correct pi/2 angle.
angles = atan(tangents);
% reduced_angles = abs(pi/3 - mod(angles,pi/3));

if ~SP_uniqueness_flag
    reduced_angles = pi/6 - abs(mod(angles,pi/3) - pi/6);
    amplitudes = sum(displacement_field.^2,3).^0.5;
    
    % see the piece of paper on the kitchen table for the math
    Values = amplitudes.*(2*cos(reduced_angles)/hlc);
    Hues = reduced_angles/(pi/6)*AB_hsv(1) + (1-reduced_angles/(pi/6))*SP_hsv(1);
    Saturations = reduced_angles/(pi/6)*AB_hsv(2) + (1-reduced_angles/(pi/6))*SP_hsv(2);
    % Saturations = Saturations.*Values;
    H_weights = [];
else
    % Value calculation is the same as before.
    reduced_angles = pi/6 - abs(mod(angles,pi/3) - pi/6);
    amplitudes = sum(displacement_field.^2,3).^0.5;
    % see the piece of paper on the kitchen table for the math
    Values = amplitudes.*(2*cos(reduced_angles)/hlc);
    
    
    full_angles = atan2(displacement_field(:,:,2),displacement_field(:,:,1)); 
    in1 = full_angles >= 0 & full_angles < pi/6;
    in2 = full_angles >= pi/6 & full_angles < pi/3;
    in3 = full_angles >= pi/3 & full_angles < pi/2;
    in4 = full_angles >= pi/2 & full_angles < 2*pi/3;
    in5 = full_angles >= 2*pi/3 & full_angles < 5*pi/6;
    in6 = full_angles >= 5*pi/6 & full_angles <= pi;
    H1 = full_angles/(pi/6)*AB_hsv(1,1) + (1-full_angles/(pi/6))*SP_hsv(1,1);
    H2 = (1-(full_angles-pi/6)/(pi/6))*AB_hsv(2,1) + ((full_angles-pi/6)/(pi/6))*SP_hsv(2,1);
    H3 = ((full_angles-pi/3)/(pi/6))*AB_hsv(2,1) + (1-(full_angles-pi/3)/(pi/6))*SP_hsv(2,1);
    H4 = (1-(full_angles-pi/2)/(pi/6))*AB_hsv(3,1) + ((full_angles-pi/2)/(pi/6))*SP_hsv(3,1);
    H5 = ((full_angles-2*pi/3)/(pi/6))*AB_hsv(3,1) + (1-(full_angles-2*pi/3)/(pi/6))*SP_hsv(3,1);
    H6 = (1-(full_angles-5*pi/6)/(pi/6))*AB_hsv(1,1) + ((full_angles-5*pi/6)/(pi/6))*SP_hsv(1,1);
    S1 = full_angles/(pi/6)*AB_hsv(1,2) + (1-full_angles/(pi/6))*SP_hsv(1,2);
    S2 = (1-(full_angles-pi/6)/(pi/6))*AB_hsv(2,2) + ((full_angles-pi/6)/(pi/6))*SP_hsv(2,2);
    S3 = ((full_angles-pi/3)/(pi/6))*AB_hsv(2,2) + (1-(full_angles-pi/3)/(pi/6))*SP_hsv(2,2);
    S4 = (1-(full_angles-pi/2)/(pi/6))*AB_hsv(3,2) + ((full_angles-pi/2)/(pi/6))*SP_hsv(3,2);
    S5 = ((full_angles-2*pi/3)/(pi/6))*AB_hsv(3,2) + (1-(full_angles-2*pi/3)/(pi/6))*SP_hsv(3,2);
    S6 = (1-(full_angles-5*pi/6)/(pi/6))*AB_hsv(1,2) + ((full_angles-5*pi/6)/(pi/6))*SP_hsv(1,2);
    
    Hues = H1.*in1 + H2.*in2 + H3.*in3 + H4.*in4 + H5.*in5 + H6.*in6;
    Saturations = S1.*in1 + S2.*in2 + S3.*in3 + S4.*in4 + S5.*in5 + S6.*in6;
    H_weights = {H1,H2,H3,H4,H5,H6};
end

Values = Values.^gamma;

% If any are out of the bounds of the hexagon, set color to white
mask = Values > 1;
Values(mask) = 1;
Hues(mask) = 0;
Saturations(mask) = 0;

Hues = mod(Hues,1);


HSV_color_stack = cat(3,Hues,Saturations,Values);
RGB_color_stack = hsv2rgb(HSV_color_stack);


end

