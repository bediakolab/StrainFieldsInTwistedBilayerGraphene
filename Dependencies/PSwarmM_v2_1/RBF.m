function [Success,Problem,Population]=RBF(Problem, Population, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The RFB function
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

% aivaz@dps.uminho.pt 16-06-2009
% lnv@mat.uc.pt       16-06-2009

% No successful step at begining
Success=0;


% Reset data from previous model
%RBFData=[];
%setappdata(Problem.RBFDataHandle,'RBFData',RBFData);

% Check for enough available point in the cache
if((Problem.CacheData.Counter<Problem.Variables+2 && ~Problem.RBFPoints) ||...
        (sum(Problem.CacheData.step)<Problem.Variables+2 && Problem.RBFPoints))
    return;
end

% We have enough points to build the model
Problem.Stats.RBF.nModels=Problem.Stats.RBF.nModels+1;

% Get Y and fY values for model build
nmax=(Problem.Variables+1)*(Problem.Variables+2)/2;

if(Problem.CacheData.Counter<=nmax)
    np=Problem.CacheData.Counter;
    Y = [Problem.CacheData.Point(1:np).x];
    fY = [Problem.CacheData.Point(1:np).fx];
else
    np=nmax;
    Far=ceil(0.2*nmax); % 20 percent far far away :)
    Near=nmax-Far;      % remaining 80 percent nearby
    start=max(1,Problem.CacheData.Counter-50*(Problem.Variables+1));
    Y = [Problem.CacheData.Point(start:start+Far-1).x...
        Problem.CacheData.Point(Problem.CacheData.Counter-Near+1:Problem.CacheData.Counter).x];
    fY = [Problem.CacheData.Point(start:start+Far-1).fx...
        Problem.CacheData.Point(Problem.CacheData.Counter-Near+1:Problem.CacheData.Counter).fx];
end

%meanfY=mean(fY);
%fY(fY>meanfY)=meanfY;

% Build A
% Preallocating memory
A=zeros(np,np);
% A is symmetric
for i=1:np
    for j=i+1:np
        A(i,j) = norm(Y(:,i)-Y(:,j))^3;
    end
    A(i+1:np,i)=A(i,i+1:np);
end


A = [ A [ ones(np,1) Y'] ;
    [ ones(np,1) Y']' zeros(Problem.Variables+1,Problem.Variables+1) ];
rhs = [ fY' ; zeros(Problem.Variables+1,1) ];

% solve the linear system
if(condest(A)>eps && condest(A)<1/eps)
    lambda = A\rhs;
else
    [U,S,V]=svd(A);
    s = diag(S);
    e = zeros(length(s),1);
    ind = s/max(abs(s)) >= eps;
    e(ind) = 1./s(ind);
    lambda = V*diag(e)*U'*rhs;
end

if(~isfinite(lambda))
    return;
end

if(Problem.PollSort==1)
    % Save data for future use
    RBFData.Y=Y;
    RBFData.lambda=lambda;
    RBFData.n=Problem.Variables;
    RBFData.np=np;
    setappdata(Problem.RBFDataHandle,'RBFData',RBFData);
end

switch Problem.TRType
    case 0
        %Thrust region
        if(Problem.Stats.RBF.nModelsPreviousSuccess>0)
            TrustSize=Population.Delta*2;
        else
            TrustSize=Population.Delta;
        end
    case 1
        TrustSize=Population.RBFDelta;
    otherwise
        error('Unknown Trust region strategy for RBF model');        
end


% Initial guess for DCA
%x = zeros(Problem.Variables,1);
x = Population.y(Population.Leader,:)';

% Trust region bounds
lb=max(Problem.LB,x-TrustSize*ones(Problem.Variables,1));
ub=min(Problem.UB,x+TrustSize*ones(Problem.Variables,1));

[f,gf] = RBFObjFun(x,Y,lambda,Problem.Variables,np);
finitial=f;

if(norm(gf)<1d-6)
    % Initial guess is stationary point
    % Do a random preturbation
    x = x + (ub-lb).*(rand(Problem.Variables,1)-0.5);
    [f,gf] = RBFObjFun(x,Y,lambda,Problem.Variables,np);
    %    display('Randomly perturbing initial guess');
end


switch Problem.RBFAlgo
    case 1 % DCA

        % DCA2 definitions
        kmax = 3000;
        rhomax = max(abs(lambda(1:np)))*6*2*max(norm(lb),norm(ub))*np;
        rho = Problem.DCARhoFactor*rhomax;
        % End of DCA2 definitions

        relerror = 1;
        k = 0;

        while (relerror>1d-6 && k<kmax)

            stop=0;
            while(~stop)
                y = rho*x - (gf - lambda(np+2:np+Problem.Variables+1));

                % Projection over the feasible set
                yaux = (y - lambda(np+2:np+Problem.Variables+1))/rho;
                masku = ( yaux > ub ); yaux(masku) = ub(masku); % Thrust region
                maskl = ( yaux < lb ); yaux(maskl) = lb(maskl);
                % End of projection to the feasible set

                [fyaux,gfyaux] = RBFObjFun(yaux,Y,lambda,Problem.Variables,np);
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

    case 2

        options = optimset('GradObj','on','Display','off','Diagnostics','off','TolFun',1d-5,'TolX',1d-5);
        [x,f]=fmincon('RBFObjFun',x,[],[],[],[],lb,ub,[],options, Y, lambda, Problem.Variables, np);
        
    case 3
        
        error('In order to use DIRECT as the RBF model minimizer you need to obtain the DIRECT solver and to comment this error message');

        DirectProblem.f='RBFObjFun';
        DirectBounds=zeros(2,length(lb))';
        DirectBounds(:,1)=lb;
        DirectBounds(:,2)=ub;
        DirectOpts.maxevals=2000;      % maximum number of function evaluations
        DirectOpts.maxits=2000;       % disable maximum number of iterations
        DirectOpts.maxdeep=100;        % maximum number of rectangle divisions
        DirectOpts.showits=0;

        [f,x]=Direct(DirectProblem,DirectBounds,DirectOpts, Y, lambda, Problem.Variables, np);

    otherwise
        error('Unknown algorithm for RBF minimization');
end

%plotRBF(Problem.Stats.RBF.nModels, lb, ub, Y,lambda,Problem.Variables,np,x);

% Damp point if linear constraints are present
if Problem.mLinear>0
    x=DampingLinearConstraints(Problem,Population,x);
end


% Evaluate true objective function
[Problem,ObjValue] = PenaltyEval(Problem, x', 0, varargin{:});

% Did we have any progress in the model?
switch Problem.TRType
    case 0
        if(Population.fy(Population.Leader)>ObjValue)
            % Success
            Success=1;
            Population.y(Population.Leader,:)=x';
            %        display('Progress');
            %        Population.fy(Population.Leader)-ObjValue
            Population.fy(Population.Leader)=ObjValue;

            % We have success in the model
            Problem.Stats.RBF.Success=Problem.Stats.RBF.Success+1;
        end
    case 1
        if(finitial-f>0)

            TRRho=(Population.fy(Population.Leader)-ObjValue)/(finitial-f);

            if(TRRho>=Problem.RBFEta1)
                Population.RBFDelta=min(Problem.RBFGamma1*Population.RBFDelta,Problem.RBFDeltaMax);
            elseif (TRRho<Problem.RBFEta0)
                Population.RBFDelta=Problem.RBFGamma0*Population.RBFDelta;
            end

            if(TRRho>Problem.RBFEta0)
                % Success
                Success=1;
                Population.y(Population.Leader,:)=x';
                %        display('Progress');
                %        Population.fy(Population.Leader)-ObjValue
                Population.fy(Population.Leader)=ObjValue;

                % We have success in the model
                Problem.Stats.RBF.Success=Problem.Stats.RBF.Success+1;
            end
        end
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
