function CacheData=InitCache(Problem)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Subrotine InitCache
%
%  This subroutine initializes the cache structure
%
%
%  Input:
%    Problem - The problem structure
%
%  Output:
%    CacheData - Cache data structure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(~isfield(Problem, 'LoadCache'))
    CacheData=[];
    return;
end

% Assume load cache is not correct
correct=0;
% A Cache file exists and will be used
if Problem.LoadCache && (exist(Problem.CacheFile,'file') || exist([char(Problem.CacheFile),'.mat'],'file'))
    if(Problem.IPrint>0)
        fprintf('Loading cache from file %s\n', char(Problem.CacheFile));
    end
    load(Problem.CacheFile,'CacheData');
    
    % check for cache correctness
    try
        correct=1;
        if ~isfield(CacheData,'xnorm')
            correct=0;
        end
        if ~isfield(CacheData,'step')
            correct=0;
        end
        if ~isfield(CacheData,'Counter')
            correct=0;
        end
        if ~isfield(CacheData,'Point')
            correct=0;
        else
            if ~isfield(CacheData.Point,'x')
                correct=0;
            end
            if ~isfield(CacheData.Point,'fx')
                correct=0;
            end
        end
    catch
        correct=0;
    end
end

% number of hits
CacheData.hits=0;
% Tolerance for compare
CacheData.Tolerance=Problem.Tolerance/10;

if(~correct)
    % number of points in the cache
    CacheData.Counter=0;
    % allocate memory in chunks
    Idx=1:Problem.CacheChunks;
    [CacheData.Point(Idx).x]=deal(zeros(Problem.Variables,1));
    [CacheData.Point(Idx).fx]=deal(0);
    [CacheData.xnorm(Idx)]=deal(0);
    [CacheData.step(Idx)]=deal(0);
end
