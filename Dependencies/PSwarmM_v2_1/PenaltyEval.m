function [Problem,ObjValue] = PenaltyEval(Problem, x, step, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Subrotine PenaltyEval
%  This subrotine is called to evaluate the objective function
%
%  Input:
%    Problem - The problem structure
%    x - Real part of point to evaluate
%    varargin - Extra parameters for objective function
%
%  Output:
%    Problem - The problem structure updated
%    ObjValue - Objective function value. Returns +Inf for unfeasible
%       points
%
% Copyright (C) 2007 A. Ismael F. Vaz and L. N. Vicente.
%
% http://www.norg.uminho.pt/aivaz
% http://www.mat.uc.pt/~lnv
%
% This library is free software; you can redistribute it and/or
% modify it under the terms of the GNU Lesser General Public
% License as published by the Free Software Foundation; either
% version 2.1 of the License, or (at your option) any later version.
%
% This library is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% Lesser General Public License for more details.
%
% You should have received a copy of the GNU Lesser General Public
% License along with this library; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%aivaz@dps.uminho.pt 06/12/2006

% Bound feasibility is enforced by the projection onto the bound feasible
% domain
% LB <= x' && x'<=UB
% Be paranoic


% There is no need to check for Problem.Vectorized, since we are assuming a
% matrix with NPoints for penalty function evaluation
[NPoints,Vars]=size(x);
ObjValue(1:NPoints)=+Inf;

% calls to the barrier/penalty function
% we account for each point
Problem.Stats.PenaltyFunCounter=Problem.Stats.PenaltyFunCounter+NPoints;


% Check for indices of points that are bound feasible
FeasibleBoundsIdx=find(min(x'-repmat(Problem.LB,1,NPoints))>=0 & min(repmat(Problem.UB,1,NPoints)-x')>=0);

if Problem.mLinear==0
	FeasibleLinearIdx=FeasibleBoundsIdx;
else
    FeasibleLinearIdx=find(max(Problem.A*x(FeasibleBoundsIdx,:)'-repmat(Problem.b,1,length(FeasibleBoundsIdx)))<=0);
end

PointsEval=x(FeasibleLinearIdx,:)';
[Vars,NPoints]=size(PointsEval);

Problem.Stats.ObjFunCounter=Problem.Stats.ObjFunCounter+NPoints;

% Check for cache hits on points
if(Problem.Cache && ~isempty(PointsEval))
    for i=NPoints:-1:1
        [yes,fx,Problem.CacheData]=CheckCache(Problem.CacheData,x(FeasibleLinearIdx(i),:)');
        if(yes)
            % its in cache, therefore we return a cached obj value
            ObjValue(FeasibleLinearIdx(i))=fx;
%            if(ObjValue(FeasibleLinearIdx(i))>10^-5+feval(Problem.ObjFunction, x(FeasibleLinearIdx(i),:)', varargin{:})...
%                    || ObjValue(FeasibleLinearIdx(i))<-10^-5+feval(Problem.ObjFunction, x(FeasibleLinearIdx(i),:)', varargin{:}))
%                ObjValue(FeasibleLinearIdx(i))
%                feval(Problem.ObjFunction, x(FeasibleLinearIdx(i),:)', varargin{:})
%                error('POOOOOOOO');
%            end
            PointsEval(:,FeasibleLinearIdx(i))=[];
            FeasibleLinearIdx(i)=[];
            NPoints=NPoints-1;
        end
    end
end

% Compute objective function value for remaining points
if(~isempty(FeasibleLinearIdx))
    try
        ObjValue(FeasibleLinearIdx)=feval(Problem.ObjFunction, PointsEval, varargin{:});
        % update counter
        Problem.Stats.RealObjFunCounter=Problem.Stats.RealObjFunCounter+NPoints;
        % Insert points in cache
        if(Problem.Cache)
            for i=1:NPoints
                Problem =...
                    InsertCache(Problem,PointsEval(:,i),ObjValue(FeasibleLinearIdx(i)),step);
            end
        end
        
    catch
        error('pswarm:ObjectiveError', ...
            ['Cannot continue because user supplied objective function' ...
            ' failed with the following error:\n%s'], lasterr)
    end
end

return;


function [yes,fx,CacheData]=CheckCache(CacheData,x)
% Check if point x is in cache.

xnorm=norm(x,inf);

% exclude points far away in norm
Idx = find((abs(CacheData.xnorm(1:CacheData.Counter) - xnorm) < CacheData.Tolerance));
if isempty(Idx)
    yes=0;
    fx=+Inf;
    return;
end

% exclude points far away (norm 1 of the distance)
points = [CacheData.Point(Idx).x];
Idx = Idx(max(abs(points - repmat(x,1,length(Idx))),[],1) < CacheData.Tolerance);
yes = ~isempty(Idx);
% returning first in cache
if(yes)
    fx= CacheData.Point(Idx(1)).fx;
else
    fx= +Inf;
end

CacheData.hits=CacheData.hits+1;

return;


function Problem = InsertCache(Problem,x,fx,step)
% Insert point in cache

% Check cache size. If needed alocate a new chunck of memory
if(Problem.CacheData.Counter>=length(Problem.CacheData.xnorm))
    % allocate a new chunck
    Idx=1+length(Problem.CacheData.xnorm):Problem.CacheChunks+length(Problem.CacheData.xnorm);
    [Problem.CacheData.Point(Idx).x]=deal(zeros(Problem.Variables,1));
    [Problem.CacheData.Point(Idx).fx]=deal(0);
    [Problem.CacheData.xnorm(Idx)]=deal(0);
    [Problem.CacheData.step(Idx)]=deal(0);
end


% Insert point
Problem.CacheData.Counter=Problem.CacheData.Counter+1;

xnorm=norm(x,inf);
Problem.CacheData.Point(Problem.CacheData.Counter).x=x;
Problem.CacheData.Point(Problem.CacheData.Counter).fx=fx;
Problem.CacheData.xnorm(Problem.CacheData.Counter)=xnorm;
Problem.CacheData.step(Problem.CacheData.Counter)=step;

return;
