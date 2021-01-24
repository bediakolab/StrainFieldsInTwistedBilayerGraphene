function [ best_DSC, convergence_storage, video_results ] = extendedMultistartDiskOptimization( objfcn, lb, ub, options, trigfun_flag, video_flag )
% Nathanael Kazmierczak, 03/22/2020
%
% Function for exerting more minute control over the optimization procedure

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
% initial_guesses = [0,0;
%                    1,0;
%                    0,1;
%                    -1,0;
%                    0.7,0.7;
%                    -0.7,0.7];
               
               initial_guesses = [0,0.05;
                   1.15,0.05;
                   0,1.3;
                   -1.15,0.05;
                   0.6,1;
                   -0.6,1;
                   1.15,0.65;
                   -1.15,0.65;
                   -0.7,0.25;
                   0.7,0.25;
                   0.3,0.65;
                   -0.3,0.65];
end
               
% hexagon_lattice_constant = 2.461;
% grid_density = 0.4;
% [ initial_guesses ] = getDSCHexagonRaster(grid_density,hexagon_lattice_constant);
               
DSC_storage = zeros(size(initial_guesses));      
RMSR_storage = zeros(size(initial_guesses,1),1);
video_results = cell(1,12);
for i = 1:size(initial_guesses,1)
    thisx0 = initial_guesses(i,:);
    successful = 0;
    count = 0;
    if video_flag
        movie_storage_null = zeros(0,2);
        assignin('base','DISPLACEMENT__STORAGE',movie_storage_null);
        options.OutputFcn = @lsqnonlin_outfun;
        % This will populate DISPLACEMENT__STORAGEe
    end
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
    if video_flag
        video_results{i} = evalin('base','DISPLACEMENT__STORAGE');
    end
end

[val,idx] = min(RMSR_storage);
best_DSC = DSC_storage(idx,:);

% Sort before storage
augmat = [RMSR_storage,DSC_storage];
sortmat = sortrows(augmat);
% DSC convergence locations first as rows of a matrix, RMSR values second. 
convergence_storage = {sortmat(:,2:end),sortmat(:,1)};

if ~video_flag
    video_results = [];
end


end

