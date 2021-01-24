function [ best_DSC ] = multistartDiskOptimization( objfcn, lb, ub, options, trigfun_flag )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% 


if ~trigfun_flag || nargin < 5
initial_guesses = [0    ,   0;
                   1    ,   0;
                   -1   ,   0;
                   0    ,   1;
                   0    ,  -1;
                   0.7  ,0.7;
                   -0.7,0.7;
                   0.7,-0.7;
                   -0.7, -0.7];
else
initial_guesses = [0,0;
                   1,0;
                   0,1;
                   -1,0;
                   0.7,0.7;
                   -0.7,0.7];
end
               
% hexagon_lattice_constant = 2.461;
% grid_density = 0.4;
% [ initial_guesses ] = getDSCHexagonRaster(grid_density,hexagon_lattice_constant);
               
DSC_storage = zeros(size(initial_guesses));      
RMSR_storage = zeros(size(initial_guesses,1),1);
for i = 1:size(initial_guesses,1)
    thisx0 = initial_guesses(i,:);
    successful = 0;
    count = 0;
    while ~successful && count < 3
        try
            DSC_storage(i,:) = lsqnonlin(objfcn,thisx0,lb,ub,options);
            successful = 1;
        catch
            successful = 0;
            count = count + 1;
        end
    end
    RMSR_storage(i) = rms(objfcn(DSC_storage(i,:)));        
end

[val,idx] = min(RMSR_storage);
best_DSC = DSC_storage(idx,:);


end

