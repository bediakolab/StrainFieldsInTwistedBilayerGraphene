function [Problem, Population]=InitPatternSearch(Problem, Population, Options, DefaultOpt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Subrotine InitPatternSearch
%
%  This subroutine initializes the search directions for the poll step
%
%  User provided directions can be included in the
%       Problem.Poll.SearchDirections.UserSpecified matrix
%  The user provided search directions will be included prior to any other
%  directions.
%
%  Provides the coordinate directions
%
%  Input:
%    Problem - The problem structure
%    Options - user defined options
%    DefaultOpt - default options
%
%  Output:
%    Problem - The problem structure updated
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% aivaz@dps.uminho.pt 01/12/2006
% updated             03/06/2009

% Compute minimum delta for pattern search
% The minimum delta available should be at least twice the Tolerance
I=intersect(find(Problem.UB<Inf),find(Problem.LB>-Inf));
if ~isempty(I)
    InitialDelta=min(Problem.UB(I)-Problem.LB(I));
    InitialDelta = InitialDelta / Problem.InitialDeltaFactor;
else
    try
        InitialDelta = Options.InitialDelta;
    catch
        InitialDelta = DefaultOpt.InitialDelta;
    end
end

% A too small initial delta?
if InitialDelta < sqrt(2*Problem.Tolerance)
    InitialDelta = sqrt(sqrt(Problem.Tolerance));
end

% Set initial delta
Problem.InitialDelta=InitialDelta;


% compute variables scaling
FiniteVars=find(isfinite(Problem.LB.*Problem.UB));
FiniteConst=find(isfinite(Problem.b));

Problem.Scale=1;
if ~isempty(FiniteVars)
    Problem.Scale = max([max(abs(Problem.UB(FiniteVars)-Problem.LB(FiniteVars))), Problem.Scale]);
end
if ~isempty(FiniteConst)
    Problem.Scale = max([max(abs(Problem.b(FiniteConst))), Problem.Scale]);
end


if Problem.IPrint>=0
    disp('Initial delta for pattern search ');
    disp(Problem.InitialDelta);
end

% Delta for the pattern search.
Population.Delta=Problem.InitialDelta;


% Keep track of what direction had success in the last poll step. The Delta
% parameter is incremented only of a success occours twice in the same
% direction.
Problem.Poll.LastSuccess=[];


% Coordinate Search for Variables
% May not be used, but it is never changed in the algorithm
Problem.Poll.SearchDirections.Base=[ones(Problem.Variables,1)...
                                   -ones(Problem.Variables,1)...
                                   eye(Problem.Variables) ...
                                   -eye(Problem.Variables)];
                               
% Linear search directions
% Is computed during the algorithm. It's pointless to provide one.
Problem.Poll.SearchDirections.Linear=[];

% A set of directions specified by user
% The computed search directions are appended to this one. If you know a
% better directions order then you should provide it here.
Problem.Poll.SearchDirections.UserSpecified =[];

return;
