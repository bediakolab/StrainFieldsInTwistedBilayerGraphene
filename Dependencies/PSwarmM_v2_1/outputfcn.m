function [res]=outputfcn(Problem, Population)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  User function called each Option.IPrint iterations
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isstruct(Problem) || ~isstruct(Population)
    return;
end

if Problem.Stats.IterCounter<=0 % first call then print header
    fprintf('\n  Iter   Act     Leader     Objective    Delta');
    fprintf('\n  ----------------------------------------------\n');
end

fprintf('    %4i   %3i   %4i   %4.6e %4.6e\n', Problem.Stats.IterCounter, ...
    Population.ActiveParticles, Population.Leader,...
    Population.fy(Population.Leader), Population.Delta);

% return a negative value to stop PSwarm execution
res=1.0;

return;