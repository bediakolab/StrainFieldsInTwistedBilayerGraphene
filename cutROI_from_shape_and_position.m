function [ roi_mask_cut ] = cutROI_from_shape_and_position( DP, position, shape )
% Goal here is to thicken up the polygon kernel defined previously as the
% cut shape. Shape comes in as a mask
%
% This is conceptually a little different from the other ROI builds,
% because here we are really just masking with a shape that was developed
% for some other region --- not a true ROI. This must be done because the
% positions are only for the original DP shape.

shape_com = getCOMofMASK( logical(shape) );
target_mat = zeros(size(DP));  % will receive the cutting mask, which currently exists in some random space
% Position needs to be in (I,J) coords, which it usually is the way the
% registration works.
shift = position - round(shape_com);  % according to formula worked out on paper.

% Move elements of the shape over, and if something goes wrong, we have
% gone out of bounds of the array and don't need to take an action, just
% like with the shifting function.
for i = 1:size(shape,1)
    for j = 1:size(shape,2)
        test1 = i+shift(1) > 0 && i+shift(1) <= size(target_mat,1);
        test2 = j+shift(2) > 0 && j+shift(2) <= size(target_mat,2);
        if test1 && test2
            target_mat(i+shift(1),j+shift(2));
            target_mat(i+shift(1),j+shift(2)) = shape(i,j);
        end
% %         try
% %             
% %             disp('No caught exception');
% %         catch
% %             % do nothing
% %             disp('caught execption');
% %         end
    end
end

% Now that the mask has been built, perform the cut.
roi_mask_cut = DP;
roi_mask_cut(~logical(target_mat)) = 0;
% If we want this to be logical, we can do so later.

end
