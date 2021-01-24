function [Problem,Population]=UpdateInfo(Problem, Population)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Subrotine UpdateInfo
%
%  This subroutine is called after the search and poll step
%
%
%  Input:
%    Problem - The problem structure
%    Population- The population structure
%
%  Output:
%    Problem - The problem structure updated
%    Population- The population structure updated
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% aivaz@dps.uminho.pt 03-06-2009


switch Problem.SearchType
    case 0 % empty search step
        % nothing to do
    case 1 % Particle swarm
        [Problem,Population]=UpdateSwarmInfo(Problem, Population);
    case {2,3} % RBF or MFN
        
    otherwise
        error('Unknown search step type')
end

return;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Subrotine UpdateInfo
%
%  This subroutine is called after the search and poll step
%  Particle swarm particles update
%
%
%  Input:
%    Problem - The problem structure
%    Population- The population structure
%
%  Output:
%    Problem - The problem structure updated
%    Population- The population structure updated
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Problem,Population]=UpdateSwarmInfo(Problem, Population)
% Compute inercia.
Inercia = Problem.InerciaInitial - ...
    (Problem.InerciaInitial-Problem.InerciaFinal)*...
    Problem.Stats.IterCounter/Problem.MaxIterations;

% Update velocity and new particle positions
% For all particles
for i=1:Population.Size
    % Is active?
    if Population.Active(i)
        % Update velocity for real variables
        Population.vx(i,:)=Projection((Inercia*Population.vx(i,:)+ ...
            Problem.Cognitial*unifrnd(0,ones(1,Problem.Variables)).*...
            (Population.y(i,:)-Population.x(i,:))+...
            Problem.Social*unifrnd(0,ones(1,Problem.Variables)).*...
            (Population.y(Population.Leader,:)-Population.x(i,:)))',...
            -Problem.MaxVelocityVect,Problem.MaxVelocityVect);
        % Update particle position and check bound limits
        
        
        AlphaMax=ones(Problem.Variables,1);
        % check for simple bound
        I=find(Population.vx(i,:)<0); % we could violate lower bounds
        if ~isempty(I)
            AlphaMax(I)=...
                min(AlphaMax(I),(Problem.LB(I)-Population.x(i,I)')...
                ./Population.vx(i,I)');
        end
        I=find(Population.vx(i,:)>0); % we could violate upper bounds
        if ~isempty(I)
            AlphaMax(I)=...
                min(AlphaMax(I),(Problem.UB(I)-Population.x(i,I)')...
                ./Population.vx(i,I)');
        end
        
        AlphaMax=max(zeros(Problem.Variables,1),AlphaMax);
        
        velocity=AlphaMax.*Population.vx(i,:)';
        
        AlphaMaxLinCons=1.0;
        if Problem.mLinear>0 && Problem.LinearStepSize
            I=find(Problem.A*velocity>0);
            if ~isempty(I)
                % We are not allowing a step longer than 1
                AlphaMaxLinCons=...
                    min(AlphaMaxLinCons,min((Problem.b(I)-Problem.A(I,:)*Population.x(i,:)')...
                    ./(Problem.A(I,:)*velocity)));
            end
        end
        
        if(AlphaMaxLinCons>0)
            
            Population.x(i,:)=Projection(Population.x(i,:)'+(AlphaMaxLinCons.*velocity),...
                Problem.LB(1:Problem.Variables), ...
                Problem.UB(1:Problem.Variables));
            
            if Problem.mLinear>0 && Problem.LinearStepSize>1 % truncate velocity
                Population.vx(i,:)=(AlphaMaxLinCons.*velocity)';
            end
        end
    end
end

% To compute population norm. Start with Leader.
MaxVelocity=norm(Population.x(Population.Leader));

% Reset number of active particles.
Population.ActiveParticles=0;
for i=1:Population.Size
    % Check if particle is active and we do not want to remove the
    % leader.
    if Population.Active(i) && Population.Leader~=i
        % compute particle velocity norm (for the stopping criteria.
        VelocityNorm=norm(Population.vx(i,:));
        Distance=norm(Population.x(i,:)-Population.x(Population.Leader,:));
        if Distance<Population.Delta && VelocityNorm<Population.Delta
            % Is neighbour
            Population.Active(i)=false;
        else
            MaxVelocity=MaxVelocity+VelocityNorm;
        end
    end
    % Account for active particles.
    if Population.Active(i)
        Population.ActiveParticles=Population.ActiveParticles+1;
    end
end

return;