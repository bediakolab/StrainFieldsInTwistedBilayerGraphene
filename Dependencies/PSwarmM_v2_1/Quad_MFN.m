function [Success,Problem,Population]=Quad_MFN(Problem, Population, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The Quadratic Minimum Forbenious Norm model
%
%
%
%  Input:
%    Problem - The problem structure
%    Population- The population structure
%
%  Output:
%    Problem - The problem structure updated
%    Population- The population structure updated
%    Success - if the search step led to success
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% aivaz@dps.uminho.pt 05-09-2009
% Part of the code is imported from the SID-PSM solver

% No successful step at begining
Success=0;


% Reset data from previous model
%RBFData=[];
%setappdata(Problem.RBFDataHandle,'RBFData',RBFData);

% Check for enough available point in the cache
if(Problem.CacheData.Counter<Problem.Variables+2)
    % more than n+2 points are requested to be in the cache
    return;
end


[quad,H,g] = quad_Frob(Population.y(Population.Leader,:)',...
    [Problem.CacheData.Point(max(1,Problem.CacheData.Counter-50*(Problem.Variables+1)):Problem.CacheData.Counter).x],...
    [Problem.CacheData.Point(max(1,Problem.CacheData.Counter-50*(Problem.Variables+1)):Problem.CacheData.Counter).fx],1,1,0);


if quad~=1
    return;
end

% We have enough points to build the model
Problem.Stats.MFN.nModels=Problem.Stats.MFN.nModels+1;

if(Problem.PollSort==1)
    % Save data for future use
    MFNData.H=H;
    MFNData.g=g;
    MFNData.n=Problem.Variables;
    MFNData.np=np;
    setappdata(Problem.MFNDataHandle,'MFNData',MFNData);
end

switch Problem.TRType
    case 0
        %Thrust region
        if(Problem.Stats.MFN.nModelsPreviousSuccess>0)
            TrustSize=Population.Delta*2;
        else
            TrustSize=Population.Delta;
        end
    case 1
        TrustSize=Population.MFNDelta;
    otherwise
        error('Unknown Trust region strategy for MFN model');
end

x = zeros(Problem.Variables,1);

switch Problem.MFNAlgo
    case 1 % DCA

        f=0;
        gf=zeros(Problem.Variables,1);
        
        lb=-TrustSize*ones(Problem.Variables,1);
        ub=TrustSize*ones(Problem.Variables,1);
        
        % DCA2 definitions
        kmax = 3000;
        rhomax = max(eig(H));
        rho = Problem.DCARhoFactor*rhomax;
        % End of DCA2 definitions

        relerror = 1;
        k = 0;

        while (relerror>1d-6 && k<kmax)

            stop=0;
            while(~stop)
                y = rho*x - (H*x+g);

                % Projection over the feasible set
                yaux = y/rho;
                masku = ( yaux > ub ); yaux(masku) = ub(masku); % Thrust region
                maskl = ( yaux < lb ); yaux(maskl) = lb(maskl);
                % End of projection to the feasible set

                [fyaux,gfyaux] = QuadObjFun(yaux,H,g);
                if(fyaux<f)
                    stop=1;
                else
                    if(rho>=rhomax)
                        stop=1;
                        % No progress and rho is too large. Should we abort
                        % or let the algorithm run with a different point
                        yaux=x;
                        fyaux=f;
                        gfyaux=gf;
                    else
                        rho=2*rho;
                    end
                end
            end

            relerror = norm(yaux-x);

            x = yaux;
            f=fyaux;
            gf=gfyaux;

            k = k+1;

        end

        dir=x;

    case 2

        % Trust region bounds
        lb=-TrustSize*ones(Problem.Variables,1);
        ub=TrustSize*ones(Problem.Variables,1);
        
        options = optimset('GradObj','on','Display','off','Diagnostics','off','TolFun',1d-6,'TolX',1d-6);
        dir=fmincon('QuadObjFun',x,[],[],[],[],lb,ub,[],options, H, g);
        
    case 3
        
        error('In order to use DIRECT as the quadratic model minimizer you need to obtain the DIRECT solver and to comment this error message');
        
        % Trust region bounds
        lb=-TrustSize*ones(Problem.Variables,1);
        ub=TrustSize*ones(Problem.Variables,1);

        DirectProblem.f='QuadObjFun';
        DirectBounds=zeros(2,length(lb))';
        DirectBounds(:,1)=lb;
        DirectBounds(:,2)=ub;
        DirectOpts.maxevals=2000;      % maximum number of function evaluations
        DirectOpts.maxits=2000;       % disable maximum number of iterations
        DirectOpts.maxdeep=100;        % maximum number of rectangle divisions
        DirectOpts.showits=0;

        [f,dir]=Direct(DirectProblem,DirectBounds,DirectOpts, H, g);

    case 4

        error('In order to use trust as the quadratic model minimizer you need to obtain the trust solver available in the optimization toolbox and to comment this error message');
        
        dir=trust(g,H,TrustSize); %Population.Delta * 2 * sqrt(2)

    otherwise
        error('Unknown algorithm for MFN minimization');
end

x=Population.y(Population.Leader,:)'+dir;

%plotMFN(Problem.Stats.MFN.nModels, lb, ub, H,g,Problem.Variables,np,x);

% Damp point if linear constraints are present
if Problem.mLinear>0
    x=DampingLinearConstraints(Problem,Population,x);
end


% Evaluate true objective function
[Problem,ObjValue] = PenaltyEval(Problem, x', 0, varargin{:});

% Did we have any progress in the model?
if(Population.fy(Population.Leader)>ObjValue)
    % Success
    Success=1;
    Population.y(Population.Leader,:)=x';
    %        display('Progress');
    %        Population.fy(Population.Leader)-ObjValue
    Population.fy(Population.Leader)=ObjValue;

    % We have success in the model
    Problem.Stats.MFN.Success=Problem.Stats.MFN.Success+1;
end

end



function [x]=DampingLinearConstraints(Problem,Population,x)

% x is bound feasible

dx=x-Population.y(1,:)';

AlphaMaxLinCons=1.0;

I=find(Problem.A*dx>0);
if ~isempty(I)
    % We are not allowing a step longer than 1
    AlphaMaxLinCons=...
        min(AlphaMaxLinCons,min((Problem.b(I)-Problem.A(I,:)*Population.y(1,:)')...
        ./(Problem.A(I,:)*dx)));
end

if(AlphaMaxLinCons>0)
    x=Projection(Population.y(1,:)'+(AlphaMaxLinCons.*dx),...
        Problem.LB(1:Problem.Variables), ...
        Problem.UB(1:Problem.Variables));
else
    x=Population.y(1,:)';
end
end
