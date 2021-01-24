function [ output_args ] = setFacetTwistAngles_refactored( input_args )
% Function removed from FourDSTEM_Analysis_Engine for the external
% triangulation_uncertainty script.
% 06/26/2020


a = hexagon_lattice_constant;
for i = 1:ntriangles
    these_vertices = con(i,:);
    v1 = v(these_vertices(1),:);
    v2 = v(these_vertices(2),:);
    v3 = v(these_vertices(3),:);
    % now that we have the coordinates, get all side lengths.
    l1 = sum((v1-v2).^2).^0.5*obj.scan_stepsize;
    l2 = sum((v1-v3).^2).^0.5*obj.scan_stepsize;
    l3 = sum((v2-v3).^2).^0.5*obj.scan_stepsize;
    obj.triangle_sidelengths(i,:) = [l1,l2,l3];
    meanl = mean([l1,l2,l3]);
    % NPK note 05/27/2020: the first formula I used below
    % implicitly employed the small angle approximation. The
    % following formula removes that approximation. However,
    % the small angle approximation is better than one part in
    % ten thousand for our datasets, so it gives the same
    % numbers as before.
    %                 obj.triangle_moire_angle(i) = rad2deg((a/10)/meanl);
    ause = a/10;
    radangle = 2*asin(ause/(2*meanl));
    obj.triangle_moire_angle(i) = rad2deg(radangle);
end

end

