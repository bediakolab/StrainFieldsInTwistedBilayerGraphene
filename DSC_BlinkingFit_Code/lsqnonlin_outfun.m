function [ stop ] = lsqnonlin_outfun( x, optimVals, state )
% NPK 03/26/2020
% For making multistart movies

% NPK 03/26/2020, for the purposes of visualizing the multistart
% convergence algorithm

    % DISPLACEMENT__STORAGE is the name of the base variable
    assignin('base','assigned',x);
    evalin('base','DISPLACEMENT__STORAGE(end+1,:) = assigned;');
stop = 0;


end

