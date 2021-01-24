function [Problem,Population]=InitPopulation(Problem, Population, ...
    InitialPopulation, PopSize, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Subrotine InitPopulation
%
%  This subroutine initializes population
%
%
%  Input:
%    Problem - The problem structure
%    InitPopulation - Initial guess population
%    PopSize - Population size
%    varargin - additional arguments for the objective function
%
%  Output:
%    Problem - The problem structure updated
%    Population - The population structure updated
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


switch Problem.SearchType
    case {0,2,3} % Empty and RBF search step
        if(PopSize>1)
            fprintf('With a null/RBF search step the population size must be one\n');
            fprintf('Forcing to be one...\n');
            PopSize=1;
        end
        JustOne=1; % Only one initial guess
        if(length(InitialPopulation)>1)
            fprintf('With a null/RBF search step just one initial guess can be considered\n');
            fprintf('Using the first feasible one provided...\n');
        end
        % Just use the same code as particle swarm
        [Problem, Population]=InitPopSwarm(Problem, Population,...
            InitialPopulation, PopSize, JustOne, varargin);
        [Problem,Population.fy(Population.Leader)]=...
            PenaltyEval(Problem, Population.y(Population.Leader,:), 1, varargin{:});
    case 1 % Particle swarm search step
        JustOne=0; % Consider all the initial guesses provided
        [Problem, Population]=InitPopSwarm(Problem, Population,...
            InitialPopulation, PopSize, JustOne, varargin);
        
    otherwise
        error('Unknown search step type')
end

return;



function [Problem, Population]=InitPopSwarm(Problem, Population, ...
    InitialPopulation, Size, JustOne, varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% subrotine InitPopulation
%    Randomly initialize the population
%    Include initial guesses, if provided by the user
%
% Input:
%   Problem - problem data
%   InitialPopulation - Inicial population provided by user
%   Size - Requested size for the population
%
% Output:
%   Problem - problem data (problem data may be update)
%   Population - population data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% aivaz@dps.uminho.pt 03-06-2009

% Maximum velocity among all particles velocity
Population.MaxVelocity=+Inf;

% Leader is the first one
Population.Leader=1;

% Number of particles in initial population that were accepted
nP=0;

if Problem.IPrint>=0
    disp('Generating initial population');
end

% Check if user provides a valid initial population
if ~isempty(InitialPopulation) && ~isstruct(InitialPopulation)
    error('pswarm:InitPopulation:InitialPopulation', 'Initial population must be defined in a structure.');
else
    % Copy the initial population for the population and initialize them
    % This is only done if the particle, after projection, is linear
    % feasible
    for i=1:length(InitialPopulation)
        % Particle position.
        Population.x(nP+1,:)=Projection(InitialPopulation(i).x,...
                Problem.LB(1:Problem.Variables), ...
                Problem.UB(1:Problem.Variables));
        if Problem.mLinear==0 || max(Problem.A*Population.x(nP+1,:)'-Problem.b)<=eps
            %particle is linear feasible
            nP=nP+1;
            
            % Best particle position.
            Population.y(nP,:)=Population.x(nP,:);
            % Particle velocities.
            Population.vx(nP,:)=zeros(1,Problem.Variables);
            % Particle is active at begining
            Population.Active(nP)=true;
            Population.fy(i)=+Inf;
            %            [Problem,Population.fy(nP)]=+Inf;
            %                PenaltyEval(Problem, Population.x(nP,:), varargin{:});
            
            % if we request just one initial guess then return
            if(JustOne)
                Population.ActiveParticles=1;
                Population.Leader=1;
                Population.Size=1;
                return;
            end
        end
    end
end

if (nP<length(InitialPopulation))
    warning('Only %d particles were linear feasible after projection', nP);
end

% Check for size
if nP>Size
    % User provided an initial feasible population greater than the population
    % size
    if Problem.IPrint>=0
        fprintf('Initial feasible population is greater than population size.\n');
        fprintf('Increasing population size from ');
        fprintf('%i',Size);
        fprintf(' to ');
        fprintf('%i',nP);
    end
    % Population size is increased to fit the number of initial guesses
    Population.Size=nP;
    Population.ActiveParticles=Population.Size;
    % no need to proceed
    return;
else
    % Otherwise just accept the proposed population size
    Population.Size=Size;
end

% Population size can no longer increase
Population.ActiveParticles=Population.Size;

% No need to proceed if we have all the initial guesses 
if(nP>=Population.Size)
    return;
end

% Ramdomly generate the remaining population
% They must be linear feasible

Ellipsoid=0;

if Problem.mLinear > 0
    if Problem.IPrint>=0
        disp('Generating Ellipsoid');
    end
    
    Ellipsoid=1; % we are computing an ellipsoid
    
    [A,b]=RealBounds(Problem);
    
    % We have linear constraints, so take care of the polytope
    maxiter = 80; tol1 = 1.e-8; tol2 = 1.e-6;

    if (nP<=0) || max(Problem.A*Population.x(1,:)'-Problem.b)>=eps ...
            || max(Population.x(1,:)'-Problem.LB)>=eps ...
            || max(Problem.UB-Population.x(1,:)')>=eps
        % We do not have a linear feasible approximation
        % Get one from presolve
        [msg,x0] = mve_presolve(A, b, maxiter,...
            zeros(Problem.Variables,1), tol1, 0); % initial guess with zeros
        if msg(1) ~= 's'
            % error on presolver.
            if Problem.IPrint>=0
                disp('MVE presolver error');
                disp(msg);
                disp('Trying to recover by adding fictitious bounds');
            end
            [A,b]=FicticiousBounds(Problem,A,b);
            [msg,x0] = mve_presolve(A, b, maxiter,...
                zeros(Problem.Variables,1), tol1, 1);
            if msg(1) ~= 's'
                if Problem.IPrint>=0
                    disp('MVE presolver error');
                    disp(msg);
                end
                
                trial=0;
                [lb,ub]=GetBounds(Problem); % to randomly generate initial guess into bound constraints
                while msg(1)~='s' && trial < Problem.MVETrials
                    [msg,x0] = mve_presolve(A, b, maxiter, unifrnd(lb,ub),...
                        tol1, 0);
                    if msg(1) ~= 's' && Problem.IPrint>=0
                        disp('MVE presolver error');
                        disp(msg);
                    end
                    trial=trial+1;
                end
                if msg(1) ~= 's'
                    if Problem.IPrint>=0
                        disp('MVE presolver error');
                        disp(msg);
                        disp('We are leaving the ellipsoid strategy');
                    end
                    Ellipsoid=0;
                end
            end
        end
    else
        x0=Population.x(1,:)'; % First as initial guess
    end

    if(Ellipsoid==1)
        % Get the ellipsoid
        % use again real bounds
        [A,b]=RealBounds(Problem);
        
        [x,E2,msg] = mve_solver(A, b, x0, maxiter, tol2);
        
        if ~msg
            if Problem.IPrint>=0
                disp('MVE solver error');
                disp('Trying to recover by adding fictitious bounds');
            end
            
            [A,b]=FicticiousBounds(Problem,A,b);

            [x,E2,msg] = mve_solver(A, b, x0, maxiter, tol2);

            %if ~msg
            %    disp('MVE solver error');
            %    disp('We are leaving the ellipsoid strategy');
            %    Ellipsoid=0;
            %end
        end
       
        %if(msg)
            [E,out] = chol(E2);

            if(out>0)
                if Problem.IPrint>=0
                    disp('Unable obtain ellipsoid');
                end
                Ellipsoid=0;
            else
                E = E';
            end

            % x - elipsoid center
            % E - Elipsoid matrix

            %hold off;

            %draw_ellipse([Problem.A; eye(Problem.Variables);...
            %    -eye(Problem.Variables)],...
            %    [Problem.b; Problem.UB; -Problem.LB],x0,x,E,[]);
        %else
        %    disp('We are leaving the ellipsoid strategy');
        %    Ellipsoid=0;
        %end
    end
end % End of ellipsoid computation

if ~Ellipsoid % get real bounds
    [lb,ub]=GetBounds(Problem);
end
        
for i=nP+1:Population.Size
    % Particle positions.
    if Ellipsoid==0
        if Problem.mLinear>0
            Population.x(i,:)=GetFeasiblePoint(Problem,lb,ub);
        else
            Population.x(i,:)=unifrnd(lb,ub);
        end
    else
        tmppt=unifrnd(-ones(Problem.Variables,1), ones(Problem.Variables,1));
        tmppt=tmppt/sqrt(tmppt'*tmppt);
        Population.x(i,:) = Projection(x +...
            (unifrnd(0, 1))^(1/Problem.Variables)*E*tmppt,...
            Problem.LB(1:Problem.Variables), ...
            Problem.UB(1:Problem.Variables));
        if max(Problem.A*Population.x(i,:)'-Problem.b)>0
            error('Not feasible pop. Internal error. This should not happen.');
        end
    end
    % Best particle position.
    Population.y(i,:)=Population.x(i,:);
    % Particle velocities
    Population.vx(i,:)=zeros(1,Problem.Variables);
    % Particle active or inactive
    Population.Active(i)=true;
    Population.fy(i)=+Inf;
%    [Problem,Population.fy(i)]=...
%        PenaltyEval(Problem, Population.x(i,:), varargin{:});
end


% plot population
%PlotPopulation(Problem, Population);

% plot ellipsoid
%t=0:0.1:2*pi;
%ellipsoid=E*[cos(t);sin(t)];
%plot(x(1)+ellipsoid(1,:),x(2)+ellipsoid(2,:));
 
return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Plot the population
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function []=PlotPopulation(Problem, Population)

if Problem.Variables ~= 2
    return;
end

hold on;

for j=1:Population.Size
    if Population.Active(j)
        plot(Population.y(j,1),Population.y(j,2),'ok');
    end
end
    
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Compute real bounds
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,b]=RealBounds(Problem)

if Problem.mLinear>0
    % Construct the polytope with linear constraints
    Id=eye(Problem.Variables); % Identity matrix for bound constraints
    A=[Problem.A; Id(find(Problem.UB<Inf),:); -Id(find(Problem.LB>-Inf),:)];
    b=[Problem.b; Problem.UB(find(Problem.UB<Inf)); -Problem.LB(find(Problem.LB>-Inf))];
else
    % Construct the polytope without linear constraints
    Id=eye(Problem.Variables); % Identity matrix for bound constraints
    A=[Id(find(Problem.UB<Inf),:); -Id(find(Problem.LB>-Inf),:)];
    b=[Problem.UB(find(Problem.UB<Inf)); -Problem.LB(find(Problem.LB>-Inf))];
end
return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Compute ficticious bounds
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,b]=FicticiousBounds(Problem,RealA,Realb)

A=RealA;
b=Realb;
Id=eye(Problem.Variables);

I=intersect(find(Problem.UB<Inf),find(Problem.LB<=-Inf));
% Upper bound finite and Lower bound not finite
if ~isempty(I)
    A=[A; -Id(I,:)];
    % Fictitious bound is three times the distance to origin
    o_I=ones(length(I),1);
    b=[b; -min(-100*o_I,Problem.UB(I)-3*abs(Problem.UB(I)))];
end
I=intersect(find(Problem.LB>-Inf), find(Problem.UB>=Inf));
% Lower bound finite and Upper bound not finite
if ~isempty(I)
    A=[A; Id(I,:)];
    o_I=ones(length(I),1);
    % Fictitious bound is three times the distance to origin
    b=[b; max(100*o_I,Problem.LB(I)+3*abs(Problem.LB(I)))];
end
I=intersect(find(Problem.UB>=Inf),find(Problem.LB<=-Inf)); % Both limites not finite
if ~isempty(I)
    A=[A; -Id(I,:); Id(I,:)];
    o_I=ones(length(I),1);
    b=[b; -o_I*min(-100, -10*min(Problem.LB>-Inf)); o_I*max(100, 10*max(Problem.UB<Inf))];
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Compute bounds (add ficticious if necessary)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [lb,ub]=GetBounds(Problem)

lb=Problem.LB;
ub=Problem.UB;

I=intersect(find(Problem.UB<Inf),find(Problem.LB<=-Inf));
% Upper bound finite and Lower bound not finite
if ~isempty(I)
    o_I=ones(length(I),1);
    lb(I)=min(-100*o_I,-Problem.UB(I)+3*abs(Problem.UB(I)));
end
I=intersect(find(Problem.LB>-Inf), find(Problem.UB>=Inf));
% Lower bound finite and Upper bound not finite
if ~isempty(I)
    o_I=ones(length(I),1);
    ub(I)=max(100*o_I,Problem.LB(I)+3*abs(Problem.LB(I)));
end
I=intersect(find(Problem.UB>=Inf),find(Problem.LB<=-Inf)); % Both limites not finite
if ~isempty(I)
    o_I=ones(length(I),1);
    lb(I)=o_I*min(-100, 10*min(Problem.LB>-Inf));
    ub(I)=o_I*max(100, 10*max(Problem.UB<Inf));
end

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Get Feasible point by try and error
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x]=GetFeasiblePoint(Problem,lb,ub)

for i=1:Problem.Trials
    x=unifrnd(lb,ub);
    if max(Problem.A*x-Problem.b)<=0
        % Is linear feasible
        return;
    end
end

disp('Failed to generate feasible initial guess');
error('Unable to proceed!');

return;

