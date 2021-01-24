function [BestParticle, BestParticleObj, RunData]=PSwarm(Problem, ...
    InitialPopulation, Options, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Particle swarm pattern search algorithm for global optimization.
%
% PSwarm algorithm takes advantage of the particle swarm algorithm in
% looking for a global optimum in the search phase of the pattern search
% algorithm. The pattern search enforces convergence for a local optimum.
%
%
% The pswarm MATLAB function syntax is
% [BestParticle, BestParticleObj, RunData] = PSwarm(Problem,
%                                  InitialPopulation, Options, P1, P2, ...)
%
% Input arguments:
%
%   Problem is a structure with some problem definitions
%
%     Problem.Variables        - Number of continuum real variables. If not
%                                specified then the size of the LB vector
%                                will be assigned to it.
%            .ObjFunction      - The name of the function to be called to
%                                when a objective function value is
%                                requested. This function syntax is
%                                [f]=myobj(x), where is a output argument
%                                - the objective function value - computed
%                                at (x), being x a vector with real
%                                variables.
%
%            .LB and .UB       - are column vector with the lower and upper
%                                bounds on the variables.
%
%            .A and .b         - Linear constraints Ax<=b
%
%   InitialPopulation is an array of structures.
%
%     InitialPopulation(1).x  -  Real variables initial guess to include in
%                                the initial population.
%                      (2).x  -  Second approximation for real ...
%       If the size of the structure array exceeds the population size then
%       the latter one will automatically be incremented.
%
%   Option is a structure with options to be passed to the algorithm. Call
%       pswarm('defaults') to obtain a list of the default options.
%
%   Just type opt=PSwarm('defaults') to obtain the default options for the
%       algorithm
%
%   P1, P2, ... are extra parameters passed to the objective evaluation
%       function
%
%
%
%
% Output arguments:
%
%   BestParticle              - The best particle obtained for the
%                               population (Leader)
%   BestParticleObj           - The best particle objective function
%                               obtained for the population (Leader
%                               objective function value)
%   RunData                   - Some statistics about the algorithm run
%
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

% aivaz@dps.uminho.pt 04/04/2007


% Default options
%  Cognitial - cognitial parameter in the velocity equation
%  InerciaFinalWeight - The value of the inercia parameter in the velocity
%             at last iteration.
%  InerciaInitialWeight - The value of the inercia parameter in the
%             velocity equation at first iteration.
%  MaxObj - Maximum number of objective function evaluations. Since the
%             algorithm is population based this maximum number of
%             objective function evaluation may be sligtly exceeded.
%  MaxIter - Maximum number of iterations allowed.
%  Size - Population size.
%  Social - Social parameter in the velocity equation.
%  MaxVelocityFactor - Velocity will be projected in the set
%             (UB-LB)*MaxVelocityFactor.
%  CPTolerance - Stopping criteria tolerance.
%  InitialDelta - Initial pattern search grid step.
%  DeltaIncreseFactor - Delta will be increased by this factor (on
%             successful poll steps).
%  DeltaDecreaseFactor - Delta willbe decreased by this factor (on
%             unsuccessful poll steps).
%  InitialDeltaFactor - The inicial Delta will be the min(UB-LB) times this
%             factor.
%  DegTolerance - Degenerancy tolerance.
%  LinearStepSize - Projection type on search step (particle swarm). =0 no
%               Step Size used on velocity . =1 Maximum Step Size if
%               computed. >1 Maximum Step Size if computed and velocity is
%               truncated for linear constraints.
%  EpsilonActive - Precision used to compute linear epsilon active
%               constraints.
%  TangentCone - 0 SID-PSM like version. 1 NOMAD like version. 2 SID-PSM
%               adapted (active constraints changes). 3 NullSpace version.
%               4 Simple QR factorization.
%
%  Trials - number of trials in generating initial guess.
%  IPrint - Iteration Print. If <0 no print is done. If =0 print the best
%           point found. If =n, n>0, print the best objective function
%           value each n iterations.
%  OutputFCN - Function to be called each IPrint iterations.
%  Vectorized - Call objective function in a vectorized way (=1) or point
%           by point (=0). In the vectorized version the poll step is
%           non-oportunistic, while in the not vectorized version the poll
%           step is oportunistic.
%  SearchType - search step type. 0 - no search. 1 - particle swarm. 2 -
%           RBF. 3 - Quadratic Model with Minumum Frobenius Norm.
%  Cache - 1 - cache objective function values
%  LoadCache - 1 - Load cache from file at startup
%  SaveCache - 1 - Save cache to file at end
%  CacheFile - Cache filename
%  RBFAlgo - Algorithm to be used to minimize the objective function, 1 -
%           DCA, 2 - fmincon
%  RBFPoint - 0 - Only point in the poll step to build the RBF model, 1 -
%           all cached points
%  PollSort - 0 - no sort in the poll step, 1 - sort by model
%  TRType - Trust region type - 0 - based on the grid parameter, 1 - based on ratio
%

% The DefaultOpt controls also the available options
% An option not provided in the DefaultOpt will be an invalid option
DefaultOpt = struct('Cognitial', 0.5, 'InerciaFinalWeight', 0.4,...
    'InerciaInitialWeight', 0.9, ...
    'MaxObj', 2000, 'MaxIter', 2000, 'Size', 20, 'Social', 0.5, ...
    'MaxVelocityFactor', 0.5, 'CPTolerance', 1.0e-5, 'InitialDelta', 20.0, ...
    'DeltaIncreaseFactor', 2.0, 'DeltaDecreaseFactor', 0.5, ...
    'InitialDeltaFactor', 5.0, 'IPrint', 10, 'DegTolerance', 0.001, ...
    'LinearStepSize', 1, 'EpsilonActive', 1e-1, 'TangentCone', 4,...
    'Trials', 1000, 'MVETrials', 10, 'OutputFCN', 'outputfcn', 'Vectorized', 0,...
    'SearchType', 1, 'Cache', 0, 'CacheChunks', 100, 'LoadCache', 0,...
    'SaveCache', 0, 'CacheFile', 'PSCache', 'RBFAlgo', 1, 'MFNAlgo', 1,...
    'DCARhoFactor', 5*10^-3,'RBFPoints',1,'PollSort',0,'TRType',0);

% With no arguments just print an error
if nargin==0 && nargout==0
    error('PSwarm:Arguments', 'Invalid number of arguments. Type ''help PSwarm'' to obtain help.');
end


% If just 'defaults' passed in, return the default options in X
if nargin==1 && nargout <= 1 && isequal(Problem,'defaults')
   BestParticle = DefaultOpt;
   return
end


% Check parameters consistence
if nargin < 1
  error('pswarm:AtLeastOneInput','PSwarm requests at least one imput (Problem definition).');
end

% If parameters are missing just define them as empty
if nargin < 3, Options=[];
    if nargin < 2, InitialPopulation = [];
    end
end

% Problem should be a structure
if ~isstruct(Problem)
    error('pswarm:StructProblem', 'First parameter must be a struct.');
end


% Do some sanity checkup on the user provided parameters and data
% We need an objective function
if ~isfield(Problem,'ObjFunction') || isempty(Problem.ObjFunction)
    error('pswarm:ObjMissing', 'Objective function name is missing.');
end

% and simple bound for particle swarm
if ~isfield(Problem,'LB') || ~isfield(Problem,'UB') || ~isnumeric(Problem.LB) || ~isnumeric(Problem.UB) || ...
        isempty(Problem.LB) || isempty(Problem.UB)
    error('pswarm:Bounds', 'PSwarm relies on bound constraints on all variables.');
end

% Bound arrays must have the same size
if length(Problem.LB)~=length(Problem.UB)
    error('pswarm:BoundsSize', 'Lower bound and upper bound arrays length mismatch.');
end

% Compute the number of variables
% If the user provides these number check for correctness, otherwise try to
% compute them from the size of the bound constraints
if (~isfield(Problem, 'Variables') || isempty(Problem.Variables))
    % default is to consider all variables as real
    Problem.Variables=length(Problem.LB);
end

% Check for computed number of variables
if Problem.Variables<0 || Problem.Variables>length(Problem.LB)
    error('pswarm:VariablesNumber', 'Number of variables do not agree with bound constraints');
end

% Be sure to have a column vector
ULBSize=size(Problem.LB);
if(ULBSize(1)~=Problem.Variables)
    Problem.LB=Problem.LB';
end
ULBSize=size(Problem.UB);
if(ULBSize(1)~=Problem.Variables)
    Problem.UB=Problem.UB';
end
for i=1:length(InitialPopulation)
    % Be sure to have a column vector
    XSize=size(InitialPopulation(i).x);
    if (XSize(1)~=Problem.Variables && XSize(2)~=Problem.Variables)
        error('pswarm:InitialPopulation', 'Initial Population size mismatch problem size');
    end
    if(XSize(1)~=Problem.Variables)
        InitialPopulation(i).x=InitialPopulation(i).x';
    end
end

% Do some checkup in the linear constraints if they exist
if isfield(Problem, 'A')
    if isempty(Problem.A)
        Problem.mLinear=0;
    else
        [Problem.mLinear, nVars]=size(Problem.A);
        if Problem.Variables ~= nVars;
            error('pswarm:Linear', 'Number of variables for Linear constraints do not agree.');
        end
    end

    % Check for independent term
    if ~isfield(Problem, 'b') || Problem.mLinear ~= length(Problem.b)
        error('pswarm:IndependentTerm', 'Size of Linear independent term does not agree or is not defined.');
    end
else
    % Always define linear constraints
    Problem.mLinear=0;
    Problem.A=[];
    Problem.b=[];
end




% check for invalid options (just to warn user about type errors).
if(isstruct(Options))
    InvalidOptions=setdiff(fieldnames(Options),fieldnames(DefaultOpt));
    if(~isempty(InvalidOptions))
        for i=1:length(InvalidOptions)
            fprintf('Invalid option(s) --> ');
            fprintf('%s\n',char(InvalidOptions(i)));
        end
        error('pswarm:InvalidOption', 'Invalid option(s) (see above)');
    end
end

% Initialize options. GetOption returns the user specified value, if the
% option exists. Otherwise returns the default option.

Problem.MaxIterations=GetOption('MaxIter',Options,DefaultOpt);
Problem.MaxEvals=GetOption('MaxObj',Options,DefaultOpt);
Problem.IPrint=GetOption('IPrint', Options, DefaultOpt);
Problem.OutputFCN=GetOption('OutputFCN', Options, DefaultOpt);
Problem.IncreaseDelta=GetOption('DeltaIncreaseFactor',Options,DefaultOpt);
Problem.DecreaseDelta=GetOption('DeltaDecreaseFactor',Options,DefaultOpt);
Problem.DegTolerance=GetOption('DegTolerance',Options,DefaultOpt);
Problem.Tolerance=GetOption('CPTolerance', Options, DefaultOpt);
Problem.LinearStepSize=GetOption('LinearStepSize', Options, DefaultOpt);
Problem.EpsilonActive=GetOption('EpsilonActive', Options, DefaultOpt);
Problem.TangentCone=GetOption('TangentCone', Options, DefaultOpt);
Problem.Trials=GetOption('Trials', Options, DefaultOpt);
Problem.MVETrials=GetOption('MVETrials', Options, DefaultOpt);
Problem.Vectorized=GetOption('Vectorized', Options, DefaultOpt);
Problem.SearchType=GetOption('SearchType', Options, DefaultOpt);
%Problem.InitialDelta=GetOption('InitialDelta', Options, DefaultOpt);
Problem.InitialDeltaFactor=GetOption('InitialDeltaFactor', Options, DefaultOpt);
Problem.Cache=GetOption('Cache', Options, DefaultOpt);
if(Problem.Cache)
    Problem.CacheChunks=GetOption('CacheChunks', Options, DefaultOpt);
    Problem.LoadCache=GetOption('LoadCache', Options, DefaultOpt);
    Problem.SaveCache=GetOption('SaveCache', Options, DefaultOpt);
    Problem.CacheFile=GetOption('CacheFile', Options, DefaultOpt);
end

Population.Size=GetOption('Size', Options, DefaultOpt);
Problem.PollSort=GetOption('PollSort', Options, DefaultOpt);
Problem.TRType=GetOption('TRType', Options, DefaultOpt);


% Initialize statistics counters

% Number of objective function calls. This is the number of calls to
% Problem.ObjFunction. Since particle swarm relies in a infinite penalty
% function strategy the penalty number of evaluation may be different.
Problem.Stats.ObjFunCounter=0;

% Number of penalty function calls. This is the number of calls to
% the barrier penalty function.
Problem.Stats.PenaltyFunCounter=0;

% Number of penalty function calls. This is the number of calls to
% the barrier penalty function.
Problem.Stats.RealObjFunCounter=0;

% Number of Poll steps taken
Problem.Stats.PollSteps=0;

% How many poll step had success
Problem.Stats.SuccPollSteps=0;

% How many times Degeneracy happens in poll step
Problem.Stats.Degenerate=0;


% Initialize patter search. Initialize the coordinate search directions.
[Problem, Population]=InitPatternSearch(Problem, Population, Options, DefaultOpt);





switch Problem.SearchType
    case 0 % Empty search step
        
    case 1 % Particle swarm search step
        Problem.InerciaInitial=GetOption('InerciaInitialWeight',Options,DefaultOpt);
        Problem.InerciaFinal=GetOption('InerciaFinalWeight',Options,DefaultOpt);
        Problem.Cognitial=GetOption('Cognitial',Options,DefaultOpt);
        Problem.Social=GetOption('Social',Options,DefaultOpt);
        Problem.MaxVelocityVect=(Problem.UB(1:Problem.Variables)-Problem.LB(1:Problem.Variables))*...
            GetOption('MaxVelocityFactor',Options,DefaultOpt);
        
        if(Problem.Cache && Problem.IPrint>0)
            fprintf('Cache of objective function values is enabled with');
            fprintf(' the particle swarm search step\n');
            fprintf('Performance may be reduced due to memory usage\n');
        end
    case 2 % RBF search step
        if(~Problem.Cache)
            fprintf('Cache of the objective function values must be enabled');
            fprintf(' for the RBF search step.\n Enabling it...');
            Problem.Cache=1;
        end
        
        % For DC algorithm
        Problem.Stats.RBF.nModels=0;
        Problem.Stats.RBF.nModelsPreviousSuccess=0;
        Problem.Stats.RBF.Success=0;
        Problem.RBFAlgo=GetOption('RBFAlgo', Options, DefaultOpt);
        Problem.RBFPoints=GetOption('RBFPoints', Options, DefaultOpt);
        Problem.DCARhoFactor=GetOption('DCARhoFactor', Options, DefaultOpt);
        
        % Trust Region parameters - for now not allowed to be changed by
        % user
        if(Problem.TRType==1)
            Problem.RBFEta0=Problem.Tolerance^0.5;
            Problem.RBFEta1=0.2;
            Problem.RBFGamma0=0.5;
            Problem.RBFGamma1=2;
            Problem.RBFDeltaMax=Problem.InitialDelta;
            if(Problem.RBFDeltaMax>1)
                Population.RBFDelta=Problem.RBFDeltaMax^0.5;
            else
                Population.RBFDelta=Problem.RBFDeltaMax^2;
            end
        end
        
        % Model data handler (appdata)
        if(Problem.PollSort==1)
            Problem.RBFDataHandle=0;
            
            % Reset data
            RBFData=[];
            setappdata(Problem.RBFDataHandle,'RBFData',RBFData);
        end
        
    case 3 % MFN search step
        if(~Problem.Cache)
            fprintf('Cache of the objective function values must be enabled');
            fprintf(' for the MFN search step.\n Enabling it...');
            Problem.Cache=1;
        end

        % For DC algorithm
        Problem.Stats.MFN.nModels=0;
        Problem.Stats.MFN.nModelsPreviousSuccess=0;
        Problem.Stats.MFN.Success=0;
        Problem.MFNAlgo=GetOption('MFNAlgo', Options, DefaultOpt);
        Problem.DCARhoFactor=GetOption('DCARhoFactor', Options, DefaultOpt);

        % Trust Region parameters - for now not allowed to be changed by
        % user
        if(Problem.TRType==1)
            Problem.RBFEta0=Problem.Tolerance^0.5;
            Problem.RBFEta1=0.2;
            Problem.RBFGamma0=0.5;
            Problem.RBFGamma1=2;
            Problem.RBFDeltaMax=Problem.InitialDelta;
            if(Problem.RBFDeltaMax>1)
                Population.MFNDelta=Problem.RBFDeltaMax^0.5;
            else
                Population.MFNDelta=Problem.RBFDeltaMax^2;
            end
        end

        % Model data handler (appdata)
        if(Problem.PollSort==1)
            Problem.MFNDataHandle=0;

            % Reset data
            MFNData=[];
            setappdata(Problem.MFNDataHandle,'MFNData',MFNData);
        end

    otherwise
        error('Unknown search step type')
end



% Initialize counters
% Iteration counter
Problem.Stats.IterCounter=0;
% How many consecutive iterations without success (progress in the leader
% particle).
IterUnsuccess=0;

% User was already warned of unfeasible leader?
Problem.FeasWarning=0;

% Cache initialization, before InitPopulation
Problem.CacheData=InitCache(Problem);

% Generate initial population. Include initial population provided by user
%  if available
[Problem,Population]=InitPopulation(Problem, Population, InitialPopulation, ...
    GetOption('Size', Options, DefaultOpt), varargin{:});

% Main cycle of the algorithm
% Run pattern search until no success is attained. Call for a poll step in
% the leader particle whenever no success is recorded.

% Start clock
tic;

% Stop if the maximum number of iterations or objective function
% evaluations is reached.
while(Problem.Stats.IterCounter<Problem.MaxIterations && Problem.Stats.ObjFunCounter<Problem.MaxEvals)
    
    if Problem.IPrint>0 && ...
        (Problem.Stats.IterCounter==1 || mod(Problem.Stats.IterCounter,Problem.IPrint)==0)
            res=feval(Problem.OutputFCN, Problem, Population);
            if(isnumeric(res) && res<0)
                if Problem.IPrint>=0
                    disp('Stopping due to user request');
                end
                break;
            end
    end
    
    if(CheckStoppingCriteria(Problem, Population))
        break;
    end
        
    % Increment iteration counter.    
    Problem.Stats.IterCounter=Problem.Stats.IterCounter+1;
    
    switch Problem.SearchType
        case 0 % Empty search step
            % no sucess to force a poll step
            Success=0;
        case 1 % Particle swarm search step
            [Success,Problem,Population]=ParticleSwarm(Problem,Population, varargin{:});
        case 2 % RBF search step
            [Success,Problem,Population]=RBF(Problem, Population, varargin{:});
        case 3 % MFN search step
            [Success,Problem,Population]=Quad_MFN(Problem, Population, varargin{:});
        otherwise
            error('Unknown search step type')
    end
    
    % Successful iteration?
    if ~Success
        % No success in search phase. Proceed to a poll step if possible.
        if Population.Delta >= Problem.Tolerance
            [Problem,Population]=PollStep(Problem,Population, ...
                varargin{:});
            Problem.Stats.PollSteps=Problem.Stats.PollSteps+1;
            % Reset on the number of unsuccessful iterations without a poll
            % step
            IterUnsuccess=0;
        else
            % An unsuccessful iteration without poll step
            IterUnsuccess=IterUnsuccess+1;
        end
    else
        % Success
        IterUnsuccess=0;
        % Leader changed.
        % Increase Delta. Reset on the local search.
        if Population.Delta<Problem.InitialDelta
            Population.Delta=Population.Delta*Problem.IncreaseDelta;
        end
        % check for lower bounds
        if Population.Delta<Problem.Tolerance
            Population.Delta=2*Problem.Tolerance;
        end
    end
    
    
    [Problem, Population]=UpdateInfo(Problem, Population);
        
end


% End of main cycle ...


% print final time
if Problem.IPrint>=0
    toc;
end

% Save cache file
if(Problem.Cache && Problem.SaveCache)
    CacheData=Problem.CacheData;
    save(Problem.CacheFile,'CacheData');
end

% Print if it was stopped due to the maximum of iterations or objective
% function evaluations
if Problem.Stats.IterCounter>=Problem.MaxIterations || Problem.Stats.ObjFunCounter>=Problem.MaxEvals
    if Problem.IPrint>=0
        disp('Maximum number of iterations or objective function evaluations reached');
    end
end

% Print for feasibility of the leader
if Problem.IPrint>=0 && Problem.mLinear>0
    if max(Problem.A*Population.y(Population.Leader,:)'-Problem.b)<=eps
        fprintf('Best particle is linear feasible\n');
    else
        error('Best particle is NOT linear feasible');
    end
end

% return leader position and objective function value
BestParticle=[Population.y(Population.Leader,:)]';
BestParticleObj=Population.fy(Population.Leader);
RunData=Problem.Stats;

if(isfield(Problem, 'Cache') && Problem.Cache>0)
    RunData.Cache.hits=Problem.CacheData.hits;
    RunData.Cache.size=Problem.CacheData.Counter;
end


if Problem.IPrint>=0
    % display some statistics
    disp(Problem.Stats);
end

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  subrotine GetOption
%
%  Input:
%    Option - option to get the value
%    Options - a list of options provided by user
%    DefaultOpt -  a list of default options
%
%  Output:
%    Value - The value specified by user for Option or the default
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Value]=GetOption(Option, Options, DefaultOpt)

% Check for user provided options
%if isempty(Options) || ~isstruct(Options)
    % User does not provides Options
%    Value=DefaultOpt.(Option);
%    return;
%end

% Try the option provided by user
try
    Value=Options.(Option);
catch
    Value=[];
end

% Option not provided by user
if isempty(Value)
    Value=DefaultOpt.(Option);    
end

return


