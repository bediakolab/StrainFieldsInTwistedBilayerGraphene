function [Problem,Population]=PollStep(Problem, Population, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  subroutine PollStep
%
%  Performe a poll step for the pattern search on the leader particle of
%  the particle swarm.
%
%  Input:
%    Problem - Problem structure
%    Population - The population
%
%  Output:
%    Poblem - Problem structure updated
%    Population - Population updated
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

% aivaz@dps.uminho.pt 30/03/2007

%
% Proceed with a poll step on the leader
%

% Columnwise search directions. Computes active linear
% contraints and returns the search directions (basis for tangent cone)
switch Problem.TangentCone
    case 0
        [Problem,D]=GetSearchDirectionsSIDPSM(Problem,Population,varargin{:});
    case 1
        [Problem,D]=GetSearchDirectionsNOMAD(Problem,Population,varargin{:});
    case 2
        [Problem,D]=GetSearchDirections(Problem,Population,varargin{:});
    case 3
        [Problem,D]=GetSearchDirectionsNullSpace(Problem,Population,varargin{:});
    case 4
        [Problem,D]=GetSearchDirectionsQR(Problem,Population,varargin{:});
    otherwise
        error('Unknown Tangent Cone parameter');
end

% remove equal directions or null directions
D=CleanDirections(Problem,D);

% Columnwise Leader
Leader=Population.y(Population.Leader,:)';

% sort directions
if(Problem.PollSort==1)
    D=SortDirections(Problem,Leader,D);
end


% Warn user if leader point is unfeasible. Convergence is not guaranteed
% in this case and it shouldn't happen
if  ~Problem.FeasWarning && Problem.mLinear>0 && max(Problem.A*Leader-Problem.b)>eps    
    error('Leader particle is unfeasible.');
    fprintf('\n\nPattern search can only get guaranteed convergence if Leader is feasible.\n');
    fprintf('\nSomething bad is going to happen!!\n\n\n');
    Problem.FeasWarning=1;
end

%plot(Leader(1),Leader(2),'*b');

if Problem.Vectorized
    % Since we are evaluating in a vectorized way, the algorithm is
    % non-oportunistic
    Trials=Projection(repmat(Leader,1,size(D,2))+Population.Delta*D,Problem.LB,Problem.UB);
    [Problem,ObjTrials]= PenaltyEval(Problem, Trials', 1, varargin{:});
    [Best,BestIdx]= min(ObjTrials);
    if Population.fy(Population.Leader)>Best
        % a successful poll step. Update counter.
        Problem.Stats.SuccPollSteps=Problem.Stats.SuccPollSteps+1;
        % update leader
        Population.y(Population.Leader,:)=Trials(:,BestIdx);
        Population.fy(Population.Leader)=ObjTrials(BestIdx);
        
        % Success obtained along the previous successful direction?
        if ~isempty(Problem.Poll.LastSuccess) && ...
                isequal(Problem.Poll.LastSuccess(:),D(:,BestIdx))
            % Yes. Increase Delta
            %if Population.Delta<Problem.InitialDelta
            Population.Delta=Population.Delta*Problem.IncreaseDelta;
            %end
        else
            % No. Update previous direction
            Problem.Poll.LastSuccess=D(:,BestIdx);
        end

        % Return. Successful poll step
        return;
    end
else
    % Oportunistic version
    % for all directions
    for i=1:size(D,2)
        % Compute trial point
        Trial=Projection(Leader+Population.Delta*D(:,i),Problem.LB,Problem.UB);
        % Compute objective function
        [Problem,ObjTrial]= PenaltyEval(Problem, Trial', 1, varargin{:});

        % Check for progress
        if Population.fy(Population.Leader)>ObjTrial

            %plot(Trial(1),Trial(2),'*r');

            % a successful poll step. Update counter.
            Problem.Stats.SuccPollSteps=Problem.Stats.SuccPollSteps+1;
            % update leader
            Population.y(Population.Leader,:)=Trial;
            Population.fy(Population.Leader)=ObjTrial;

            % Success obtained along the previous successful direction?
            if ~isempty(Problem.Poll.LastSuccess) && ...
                    isequal(Problem.Poll.LastSuccess(:),D(:,i))
                % Yes. Increase Delta
                %if Population.Delta<Problem.InitialDelta
                Population.Delta=Population.Delta*Problem.IncreaseDelta;
                %end
            else
                % No. Update previous direction
                Problem.Poll.LastSuccess=D(:,i);
            end

            % Return. Successful poll step
            return;
        end
    end
end

% No success in poll step. Decrease Delta.
if Population.Delta>=Problem.Tolerance
    Population.Delta=Population.Delta*Problem.DecreaseDelta;
end

Problem.Poll.LastSuccess=[];

return;


function [Problem,D]=GetSearchDirectionsSIDPSM(Problem, Population, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% compute all search directions. Coordinate and linear active
% directions - Adapted from SID-PSM
%
% Input:
%   Problem - Problem structure
%   Population - Population structure
%
% Output:
%   D - Search directions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The best population position. Columnwise
LeaderPoint=Population.y(Population.Leader,:)';

% Compute search direction with respect to near active linear constraints
if Problem.mLinear>0
    Active=Problem.A(find(Problem.A*LeaderPoint-Problem.b >= ...
        -min(Problem.EpsilonActive,10*Population.Delta)),:)';
    NActive=size(Active,2);
    
    if NActive>0 && NActive<=Problem.Variables
        % Code from Ana and Luís
        [U,S,V] = svd(Active);
        S1      = S(1:NActive,1:NActive);
        if min(diag(S1)) < Problem.DegTolerance
            Problem.Poll.SearchDirections.Linear = [];
            Problem.Stats.Degenerate=Problem.Stats.Degenerate+1;
        else
            U1  = U(1:Problem.Variables,1:NActive);
            W   = U1 * diag(1./diag(S1)) * V';
            if NActive < Problem.Variables
                U2 = U(1:Problem.Variables,NActive+1:Problem.Variables);
                d  = sum(U2',1)';
                PP=[ -U2 d -W ];
                Problem.Poll.SearchDirections.Linear=PP(1:Problem.Variables,:);
                %PP(Problem.RealVariables+1:Problem.n,:);
                %
                % Other alternatives:
                %
                %           [ U2 -d -W ];
                %
            else
                Problem.Poll.SearchDirections.Linear=-W(1:Problem.Variables,:);
            end
        end
    else
        Problem.Poll.SearchDirections.Linear=[];
    end
else
    Problem.Poll.SearchDirections.Linear=[];
end


% The direction are the concatenation of the coordinate, linear,
% and user directions
D = [Problem.Poll.SearchDirections.UserSpecified, ...
    Problem.Poll.SearchDirections.Base, ...
    Problem.Poll.SearchDirections.Linear];

return;


function [Problem,D]=GetSearchDirectionsNOMAD(Problem, Population, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% compute all search directions. Coordinate and linear active
% directions - Adapted from NOMAD
%
% Input:
%   Problem - Problem structure
%   Population - Population structure
%
% Output:
%   D - Search directions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The best population position. Columnwise
LeaderPoint=Population.y(Population.Leader,:)';

% Feasible domain in NOMAD syntax, so no change in the NOMAD code is
% requested
Omega.A = [eye(Problem.Variables), Problem.A']';
Omega.u = [Problem.UB; Problem.b];
Omega.l = [Problem.LB; Inf*ones(Problem.mLinear,1)];

try
    [B,N] = getTangentCone(LeaderPoint,Omega,Problem.EpsilonActive, ones(Problem.Variables,1));
catch
    disp('');
    disp('NOMAD version of Tangent Cone requests functions:');
    disp(' getTangentCone, scaleDirections, activeConstraints and removeRedundancy');
    disp('');
    disp('Please visit NOMAD homepage and provide these functions');
    error('NOMAD functions missing');
end



% The direction are the concatenation of the coordinate, linear,
% and user directions
D = [Problem.Poll.SearchDirections.UserSpecified N -N B -B];

return;


function [Problem,D]=GetSearchDirections(Problem, Population, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% compute all search directions. Coordinate and linear active
% directions
%
% Input:
%   Problem - Problem structure
%   Population - Population structure
%
% Output:
%   D - Search directions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The best population position. Columnwise
LeaderPoint=Population.y(Population.Leader,:)';

ActiveThreshold=min(Problem.EpsilonActive,10*Population.Delta); % should be <=1
ActiveThresholdLimit=min(0.1,ActiveThreshold^2);

% Compute search direction with respect to near active linear constraints
if Problem.mLinear>0
    while (ActiveThreshold>=ActiveThresholdLimit) % to be flexible with active constraints
        
        Active=Problem.A(find(Problem.A*LeaderPoint-Problem.b >= ...
            -ActiveThreshold*Problem.Scale),:)';
        NActive=size(Active,2);

        if NActive>0 && NActive<=Problem.Variables
            % Code from Ana and Luís
            [U,S,V] = svd(Active);
            S1      = S(1:NActive,1:NActive);
            if min(diag(S1)) < Problem.DegTolerance
                ActiveThreshold=ActiveThreshold/2;
            else
                U1  = U(1:Problem.Variables,1:NActive);
                W   = U1 * diag(1./diag(S1)) * V';
                if NActive < Problem.Variables
                    U2 = U(1:Problem.Variables,NActive+1:Problem.Variables);
                    d  = sum(U2',1)';
                    PP=[ -U2 d -W ];
                    Problem.Poll.SearchDirections.Linear=PP(1:Problem.Variables,:);
                    %PP(Problem.RealVariables+1:Problem.n,:);
                    %
                    % Other alternatives:
                    %
                    %           [ U2 -d -W ];
                    %
                else
                    Problem.Poll.SearchDirections.Linear=-W(1:Problem.Variables,:);
                end
                break; % break while
            end
        else
            Problem.Poll.SearchDirections.Linear=[]; % to be defined
            if NActive<=0
                break; % break while, since we have no active constraints
            else
                ActiveThreshold=ActiveThreshold/2; % Attempt to reduce the number of active constraints
            end
        end
    end
else
    Problem.Poll.SearchDirections.Linear=[];
end

if ActiveThreshold<ActiveThresholdLimit
    Problem.Poll.SearchDirections.Linear = [];
    Problem.Stats.Degenerate=Problem.Stats.Degenerate+1;
end


% The direction are the concatenation of the coordinate, linear,
% and user directions
D = [Problem.Poll.SearchDirections.UserSpecified,...
    Problem.Poll.SearchDirections.Base, ...
    Problem.Poll.SearchDirections.Linear];

return;



function [Problem,D]=GetSearchDirectionsQR(Problem, Population, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% compute all search directions. Coordinate and linear active
% directions
%
% Input:
%   Problem - Problem structure
%   Population - Population structure
%
% Output:
%   D - Search directions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The best population position. Columnwise
LeaderPoint=Population.y(Population.Leader,:)';

ActiveThreshold=min(Problem.EpsilonActive,10*Population.Delta); % should be <=1
ActiveThresholdLimit=min(0.1,ActiveThreshold^2);



% Compute search direction with respect to near active linear constraints
if Problem.mLinear>0
    Problem.Poll.SearchDirections.Linear=[]; % to be defined
    while (ActiveThreshold>=ActiveThresholdLimit) % to be flexible with active constraints
        
        LActive=Problem.A(find(Problem.A*LeaderPoint-Problem.b >= ...
            -ActiveThreshold*Problem.Scale),:)';
        I=eye(Problem.Variables);
        LBActive=...
            -I(find(LeaderPoint-Problem.LB<=ActiveThreshold*Problem.Scale),:);
        UBActive=...
            I(find(Problem.UB-LeaderPoint<=ActiveThreshold*Problem.Scale),:);
        
        Active=[LActive, LBActive', UBActive'];
        
        NActive=size(Active,2);

        if NActive>0 && NActive<=Problem.Variables...
                && (rank(Active) == min(size(Active)))
            [Q,R] = qr(Active,0);
            B = Q/R';
            N = I - B*Active';
            for i=Problem.Variables:-1:1
                if N(:,i)'*N(:,i)<Problem.Tolerance
                    N(:,i)=[];
                end
            end
            for i=NActive:-1:1
                if B(:,i)'*B(:,i)<Problem.Tolerance
                    B(:,i)=[];
                end
            end
                
            Problem.Poll.SearchDirections.Linear = [N, -N, B, -B];
            break;
        else
            if NActive<=0
                break; % break while, since we have no active constraints
            else
                ActiveThreshold=ActiveThreshold/2; % Attempt to reduce the number of active constraints
            end
        end
    end
else
    Problem.Poll.SearchDirections.Linear=[];
end

if ActiveThreshold<ActiveThresholdLimit
    Problem.Poll.SearchDirections.Linear = [];
    Problem.Stats.Degenerate=Problem.Stats.Degenerate+1;
end


% The direction are the concatenation of the linear,
% and user directions
if(~isempty(Problem.Poll.SearchDirections.Linear))
    D = [Problem.Poll.SearchDirections.UserSpecified,...
        Problem.Poll.SearchDirections.Linear];
else
    D = [Problem.Poll.SearchDirections.UserSpecified,...
        Problem.Poll.SearchDirections.Base];
end

return;





function [Problem,D]=GetSearchDirectionsNullSpace(Problem, Population, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% compute all search directions. Coordinate and linear active
% directions
%
% Directions are NullSpace of active constraints
%
% Input:
%   Problem - Problem structure
%   Population - Population structure
%
% Output:
%   D - Search directions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The best population position. Columnwise
LeaderPoint=Population.y(Population.Leader,:)';


% Compute search direction with respect to near active linear constraints
if Problem.mLinear>0 % We have linear constraints
    
    Active=find(Problem.A*LeaderPoint-Problem.b+...
        min(Problem.EpsilonActive,10*Population.Delta)*Problem.Scale>0);
    
    if ~isempty(Active)
        [Closest, Closesti] = min(max(Problem.A(Active,:)*LeaderPoint-Problem.b(Active)));
    else
        Closesti=0;
    end
    
    if Closesti>0 % We have an active linear constraint
        %Z = NullSpace(Problem.A(Closesti,:)');
        %Problem.Poll.SearchDirections.Linear=[Z, -Z];
        [Q,R] = qr(Problem.A(Closesti,:)',0);
        B = Q/R';
        N = eye(Problem.Variables) - B*Problem.A(Closesti,:);
        Problem.Poll.SearchDirections.Linear = [Problem.A(Closesti,:)', -Problem.A(Closesti,:)', N, -N, B, -B];
    end
else
    Problem.Poll.SearchDirections.Linear=[];
end


% The direction are the concatenation of the coordinate, linear,
% and user directions
D = [Problem.Poll.SearchDirections.UserSpecified,...
    Problem.Poll.SearchDirections.Base, ...
    Problem.Poll.SearchDirections.Linear];

return;


function [D]=CleanDirections(Problem,D)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Cleans the search directions set D by removing the null directions and
% the repeated directions
%
% Input:
%   Problem - Problem structure
%   D - Set of search directions
%
% Output:
%   D - Set of search directions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% remove any null directions 
[n,m]=size(D);
for i=m:-1:1 % decreasing order, since we may remove some directions
    if norm(D(:,i))<Problem.Tolerance
        D(:,i)=[];
    end
end

% remove repeated directions
[n,m]=size(D); %It should be faster then to update m in the previous cycle
for i=m:-1:1 %decreasing
    for j=i-1:-1:1 % decreasing
        if norm(D(:,i)-D(:,j))<Problem.Tolerance
            D(:,i)=[]; % remove direction with higher indice
            break; % break internal cycle
        end
    end
end

return;


function [D]=SortDirections(Problem,Leader,D)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Sort the search directions
%
% Input:
%   Problem - Problem structure
%   D - Set of search directions
%
% Output:
%   D - Set of search directions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch Problem.SearchType
    case 2
        % In the RBF case we are sorting by the RBF model value
        
        % get the model data, if available
        if(isappdata(Problem.RBFDataHandle,'RBFData'))
            RBFData = getappdata(Problem.RBFDataHandle,'RBFData');
            if(~isempty(RBFData))
                nDir=size(D,2);
                modelf=zeros(1,nDir);
                for i=1:nDir
                    [modelf(i),grad_f] = RBFObjFun(Leader+D(:,i), RBFData.Y, RBFData.lambda, RBFData.n, RBFData.np);
                end
                [modelfs,idx]=sort(modelf);
                D(1:Problem.Variables,1:nDir)=D(1:Problem.Variables,idx);
            end
        end
        
    case 3
        % In the MFN case we are sorting by the MFN model value

        % get the model data, if available
        if(isappdata(Problem.MFNDataHandle,'MFNData'))
            MFNData = getappdata(Problem.MFNDataHandle,'MFNData');
            if(~isempty(MFNData))
                nDir=size(D,2);
                modelf=zeros(1,nDir);
                for i=1:nDir
                    [modelf(i),grad_f] = QuadObjFun(Leader+D(:,i), MFNData.H, MFNData.g);
                end
                [modelfs,idx]=sort(modelf);
                D(1:Problem.Variables,1:nDir)=D(1:Problem.Variables,idx);
            end
        end
    otherwise
        display('Poll step with an unknowing sorting strategy');
end



return;