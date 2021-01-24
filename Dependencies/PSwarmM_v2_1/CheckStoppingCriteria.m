function [stop]=CheckStoppingCriteria(Problem, Population)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Subrotine CheckStoppingCriteria
%
%  This subroutine check for the stopping criteria
%
%
%  Input:
%    Problem - The problem structure
%    Population - The population structure
%
%  Output:
%    1 to stop
%    0 don't stop
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% aivaz@dps.uminho.pt 03-06-2009



switch Problem.SearchType
    case 0 % Empty search step
        % Stop if pattern search Delta is bellow Tolerance
        if Population.Delta<Problem.Tolerance
            if Problem.IPrint>=0
                disp('Stopping due to tolerance');
            end
            stop=1;
            return;
        end
    case 1 % Particle swarm search step
        % Stop also if particle swarm velocity is bellow Tolerance and pattern
        % search Delta is bellow Tolerance
        if Population.MaxVelocity<Problem.Tolerance && Population.Delta<Problem.Tolerance
            if Problem.IPrint>=0
                disp('Stopping due to velocity and tolerance');
            end
            stop=1;
            return;
        end
        
        % Stop if the number of active particles is equal to 1 and pattern
        % search Tolerance was attained.
        if Population.ActiveParticles <=1 && Population.Delta<Problem.Tolerance
            if Problem.IPrint>=0
                disp('Stopping due to single particle and tolerance');
            end
            stop=1;
            return;
        end
    case {2,3} % RBF or MFN search step
        % Stop if pattern search Delta is bellow Tolerance
        if Population.Delta<Problem.Tolerance
            if Problem.IPrint>=0
                disp('Stopping due to tolerance');
            end
            stop=1;
            return;
        end
    otherwise
        error('Unknown search step type')
end





stop=0;
return;