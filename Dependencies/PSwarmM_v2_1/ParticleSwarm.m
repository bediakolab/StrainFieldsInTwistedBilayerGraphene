function [Success,Problem,Population]=ParticleSwarm(Problem, Population, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Subrotine InitSwarm
%
%  This subroutine implements the particle swarm search step
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

% No success iteration at begining
Success=false;


if(Problem.Vectorized)
    % Call the penalty function for all active particles in a vectorized way
    % This is just a question of memory storage, as PenaltyEval can
    % manage a vector of points
    
    % Get all active particles in the swarm
    ActiveParticlesIdx=find(Population.Active(:)==1);
    ActiveParticles=Population.x(ActiveParticlesIdx,:);
    [Problem,ObjValue]=...
        PenaltyEval(Problem, ActiveParticles, 0, varargin{:});
    for i=1:length(ActiveParticlesIdx)
        % Was progress attained for the current particle?
        if Population.fy(ActiveParticlesIdx(i))>ObjValue(i)
            % Yes. Update best particle position
            Population.fy(ActiveParticlesIdx(i))=ObjValue(i);
            Population.y(ActiveParticlesIdx(i),:)=Population.x(ActiveParticlesIdx(i),:);
            
            % Check if new leader is available
            if Population.fy(Population.Leader)>Population.fy(ActiveParticlesIdx(i))...
                    || Population.Leader==ActiveParticlesIdx(i)
                Population.Leader=ActiveParticlesIdx(i);
                % Particle swarm iteration declared as successful
                Success=true;
                % Reset last success direction for pattern search
                Problem.Poll.LastSuccess=[];
            end
        end
    end
else
    % For all particles in the swarm
    for i=1:Population.Size
        % If particle is active
        if (Population.Active(i))
            % Compute Objective function value
            [Problem,ObjValue]=...
                PenaltyEval(Problem, Population.x(i,:), 0, varargin{:});
            % Was progress attained for the current particle?
            if Population.fy(i)>ObjValue
                % Yes. Update best particle position
                Population.fy(i)=ObjValue;
                Population.y(i,:)=Population.x(i,:);
                
                % Check if new leader is available
                if Population.fy(Population.Leader)>Population.fy(i) || Population.Leader==i
                    Population.Leader=i;
                    % Particle swarm iteration declared as successful
                    Success=true;
                    % Reset last success direction for pattern search
                    Problem.Poll.LastSuccess=[];
                end
            end
        end
    end
end

return;