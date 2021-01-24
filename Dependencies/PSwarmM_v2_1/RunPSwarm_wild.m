function [x,fx,nfo,deg,nit,npoll,spoll,nModels,RBFSuc]=RunPSwarm_wild(kp,np,mp,sp,type)
%
%
%

%try % Avoid errors in a batch run
    x=dfoxs(np,kp,10^(sp-1));
    n=size(x,1);
    Problem.LB=-Inf*ones(n,1);
    Problem.UB=+Inf*ones(n,1);
    
    % Problem definition
    Problem.ObjFunction='wild_obj';
    
    if ~isempty(x)
        InitPop(1).x=x;
    else
        InitPop=zeros(length(Problem.LB),1);
    end

    Problem.A=[];
    Problem.b=[];


    % Algorithm options
    %Options.Size=40;
    %Options.Size=1;
    %Options.MaxIter=opt;
    %Options.MaxObj=opt;
%    Options.IPrint=1;
    Options.InitialDelta=max(1,norm(InitPop(1).x,inf)); %same as SID-PSM
    Options.SearchType=0; %0- No search step, 1- Particle swarm, 2- RBF version, 3 - Quadratic model
    Options.RBFAlgo=1; % 1- DCA, 2- fmincon, 3 - Direct
    Options.MFNAlgo=1; % 1- DCA, 2- fmincon, 3 - Direct, 4 - trust
    %Options.DCARhoFactor=opt;
    %Options.PollSort=1;
    %Options.TRType=1;

    % Run the algorithm
    [x,fx,RunData]=PSwarm(Problem, InitPop, Options, kp, np, mp, sp, type);
    nfo=RunData.ObjFunCounter;
    deg=RunData.Degenerate;
    nit=RunData.IterCounter;
    npoll=RunData.PollSteps;
    spoll=RunData.SuccPollSteps;
    nModels=0; %RunData.RBF.nModels;%RunData.RBF.nModels;
    RBFSuc=0; %RunData.RBF.Success;%RunData.RBF.Success;
%catch
%end

display(x);
display(fx);

%display(RunData.MFN);