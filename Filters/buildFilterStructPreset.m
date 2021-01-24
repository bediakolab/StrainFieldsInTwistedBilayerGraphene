function [ filterstruct ] = buildFilterStructPreset( preset_id_string )
% Interfaces into buildFilterStruct.m. Most of the preset ids will be
% strings referring to the name of the dataset and the type of filtering
% being performed. For instance, "DS4 displacement #1" or "DS 26
% strain #3".
%
% Nathanael Kazmierczak, 05/19/2020

switch preset_id_string
    case 'blank filter'  % Massively permissive threshold to have no effect.
        string_cell = {'median outlier','xdisp',[3,1000]};
        filterstruct = buildFilterStruct(string_cell);
    case 'standard displacement multistart'
        string_cell = {'soft multistart','amplitude',[0.3,5]};
        filterstruct = buildFilterStruct(string_cell);
    case 'DS4 annealed displacement #1'
        stringcell = {'median outlier','xdisp',[5,0.3];'median outlier','ydisp',[5,0.3];'median outlier','xdisp',[3,0.25];'median outlier','ydisp',[3,0.25]};
        [ filterstruct ] = buildFilterStruct( stringcell );
    case 'DS4 annealed displacement #2'
        stringcell = {'median outlier','xdisp',[5,0.1];'median outlier','ydisp',[5,0.1];};
        [ filterstruct ] = buildFilterStruct( stringcell );
    case 'DS4 annealed displacement #3'
        stringcell = {'median normal','xdisp',5;'median normal','ydisp',5};
        [ filterstruct ] = buildFilterStruct( stringcell );
    case 'DS4 annealed displacement #4'
        stringcell = {'median normal','xdisp',3;'median normal','ydisp',3};
        [ filterstruct ] = buildFilterStruct( stringcell );
    case 'DS4 annealed displacement #5'
        stringcell = {'median normal','xdisp',7;'median normal','ydisp',7};
        [ filterstruct ] = buildFilterStruct( stringcell );
    case 'DS4 annealed displacement #6'
        stringcell = {'median normal','xdisp',11;'median normal','ydisp',11};
        [ filterstruct ] = buildFilterStruct( stringcell );
    case 'DS4 annealed displacement #7'
        stringcell = {'moving average','xdisp',11;'moving average','ydisp',11};
        [ filterstruct ] = buildFilterStruct( stringcell );
    case 'DS4 annealed displacement #8'
        stringcell = {'moving average','xdisp',5;'moving average','ydisp',5};
        [ filterstruct ] = buildFilterStruct( stringcell );
    case 'DS4 annealed displacement #9'
        stringcell = {'moving average','xdisp',3;'moving average','ydisp',3};
        [ filterstruct ] = buildFilterStruct( stringcell );
    case 'DS4 annealed displacement #10'
        stringcell = {'total variation','xdisp',0.01;'total variation','ydisp',0.01};
        [ filterstruct ] = buildFilterStruct( stringcell );
    case 'DS4 annealed displacement #11'
        stringcell = {'median outlier','xdisp',[5,0.3];'median outlier','ydisp',[5,0.3];'total variation','xdisp',0.01;'total variation','ydisp',0.01};
        [ filterstruct ] = buildFilterStruct( stringcell );
    case 'DS4 annealed displacement #12'
        stringcell = {'median outlier','xdisp',[5,0.3];'median outlier','ydisp',[5,0.3];'total variation','xdisp',0.05;'total variation','ydisp',0.05};
        [ filterstruct ] = buildFilterStruct( stringcell );
    case 'DS4 annealed displacement #13'
        stringcell = {'median outlier','xdisp',[5,0.3];'median outlier','ydisp',[5,0.3];'total variation','xdisp',0.5;'total variation','ydisp',0.5};
        [ filterstruct ] = buildFilterStruct( stringcell );
    case 'DS4 annealed displacement #14'
        stringcell = {'median outlier','xdisp',[5,0.2];'median outlier','ydisp',[5,0.2];'moving average','xdisp',11;'moving average','ydisp',11};
        [ filterstruct ] = buildFilterStruct( stringcell );
    case 'DS26 annealed displacement #1'
        stringcell = {'median outlier','xdisp',[5,0.2];'median outlier','ydisp',[5,0.2];'gaussian average','xdisp',5;'gaussian average','ydisp',5};
        [ filterstruct ] = buildFilterStruct( stringcell );
        % Too blobish. Probably want literal moving mean or median with
        % circular ROI. 
    case 'DS26 annealed displacement #2'
        stringcell = {'median outlier','xdisp',[5,0.2];'median outlier','ydisp',[5,0.2];'moving average circle','xdisp',11;'moving average circle','ydisp',11};
        [ filterstruct ] = buildFilterStruct( stringcell );
    case 'DS26 annealed displacement #3'
        stringcell = {'median outlier','xdisp',[5,0.2];'median outlier','ydisp',[5,0.2];'moving average circle','xdisp',21;'moving average circle','ydisp',21};
        [ filterstruct ] = buildFilterStruct( stringcell );
        
    case 'DS15 annealed displacement #1'
        stringcell = {'median outlier','xdisp',[3,0.2];'median outlier','ydisp',[3,0.2];'moving average circle','xdisp',5;'moving average circle','ydisp',5};
        [ filterstruct ] = buildFilterStruct( stringcell );
end


end

