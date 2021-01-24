function [x,fx,nfo,deg,nit,npoll,spoll,nModels,RBFSuc]=RunPswarm(prob,opt)
%
% Run a problem in the NL format from AMPL
%

%try % Avoid errors in a batch run
    [x,Problem.LB,Problem.UB,v,lc,uc,A]=matampl(prob);
    % Problem definition
    Problem.ObjFunction='ampl_obj';
    if ~isempty(x)
        InitPop(1).x=x;
    else
        InitPop=zeros(length(Problem.LB),1);
    end
    
    %InitPop=[];
%    InitPop(1).x=[10 10]';%ones(length(Problem.LB),1); % We are currently ignoring an initial guess provided by AMPL
    % Linear constraints
    I=find(lc>-Inf);
    if ~isempty(I)
        Problem.A=-A(I,:);
        Problem.b=-lc(I);
    else % define Problem.A and .b
        Problem.A=[];
        Problem.b=[];
    end
    
    I=find(uc<Inf);
    if ~isempty(I)
        Problem.A=[Problem.A; A(I,:)];
        Problem.b=[Problem.b; uc(I)];
    end
    

    % Algorithm options
    %Options.Size=40;
    %Options.Size=1;
    %Options.MaxIter=opt;
    %Options.MaxObj=opt;
    %Options.SearchType=3; %0- No search step, 1- Particle swarm, 2- RBF version, 3 - Quadratic model
    %Options.RBFAlgo=2; % 1- DCA, 2- fmincon, 3 - Direct
    %Options.MFNAlgo=4; % 1- DCA, 2- fmincon, 3 - Direct, 4 - trust
    %Options.DCARhoFactor=opt;
    %Options.PollSort=1;
    %Options.TRType=1;

    % Run the algorithm
    [x,fx,RunData]=PSwarm(Problem, InitPop, Options);
    nfo=RunData.ObjFunCounter;
    deg=RunData.Degenerate;
    nit=RunData.IterCounter;
    npoll=RunData.PollSteps;
    spoll=RunData.SuccPollSteps;
    nModels=0;%RunData.MFN.nModels;%RunData.MFN.nModels;%RunData.RBF.nModels;
    RBFSuc=0;%RunData.MFN.Success;%RunData.MFN.Success;%RunData.RBF.Success;
%catch
%end

display(x);
display(fx);

%display(RunData.RBF);