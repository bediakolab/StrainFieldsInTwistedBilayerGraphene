function [ filterstruct ] = buildFilterStruct( string_cell )
% Function for building the struct inputs to filterDisplacement.m and
% filterStrain.m. Useful for the revised filtering protocol in
% FourDSTEM_Analysis_Engine.
%
% string_cell should be an n x 3 cell array where each row specifies a
% filter to be built. The first column is .type, the second .target, the
% third .params. As a shortcut, multiple targets can be specified through
% keywords. 'four strain' will build filters with targets exx, exy, eyx, 
% and eyy. 'six strain' will also add gxy and gtorsion. 'both
% displacements' will build a filter for both x displacement and y
% displacement targets. Note that this is different than filtering the
% amplitude; the real purpose here is when doing, say, TV min on annealed
% datasets.
%
% string_cell{1}: method name ("type")
% string_cell{2}: matrix to use ("target")
% string_cell{3}: thresholds for method ("params")
%
% Nathanael Kazmierczak, 05/19/2020

nrows = size(string_cell,1);
filterstruct = struct;
filterstruct.name = [];  % empty declarations required for the assignments to work.fs
filterstruct.target = [];
filterstruct.params = [];
count = 1;
for i = 1:nrows
    % logic operates first on the target, because this will dictate if we
    % are actually making multiple filters
    target = string_cell{i,2};
    targets = cell(0,1);
    switch target
        case 'four strain'
            loop = 1;
            targets{1} = 'exx';
            targets{2} = 'exy';
            targets{3} = 'eyx';
            targets{4} = 'eyy';
        case 'six strain'
            loop = 1;
            targets{1} = 'exx';
            targets{2} = 'exy';
            targets{3} = 'eyx';
            targets{4} = 'eyy';
            targets{5} = 'gxy';
            targets{6} = 'gtorsion';
        case 'both displacements'
            loop = 1;
            targets{1} = 'xdisp';
            targets{2} = 'ydisp';
        otherwise
            loop = 0;
    end
    
    if loop 
        ntargets = numel(targets);
        for j = 1:ntargets
            this_filter = struct;
            this_filter.name = string_cell{i,1};
            this_filter.target = targets{j};
            this_filter.params = string_cell{i,3};  % Note that this should be an array, not a string
            filterstruct(count) = this_filter;
            count = count + 1;
        end
    else
        this_filter = struct;
        this_filter.name = string_cell{i,1};
        this_filter.target = target;
        this_filter.params = string_cell{i,3};
        filterstruct(count) = this_filter;
        count = count + 1;
    end
end


end

