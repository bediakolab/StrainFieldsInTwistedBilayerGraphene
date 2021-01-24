function [x,fx,nfo,deg,nit,npoll,spoll,hits,csize]=RunPswarm(prob) %,opt)
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
        InitPop=[];
    end
    
    InitPop=[]; % We are currently ignoring an initial guess provided by AMPL
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
    Options.Size=40;
    %Options.Size=1;
    %Options.MaxIter=opt;
    %Options.MaxObj=opt;

    % Run the algorithm
    [x,fx,RunData]=PSwarm(Problem, InitPop, Options);
    nfo=RunData.ObjFunCounter;
    deg=RunData.Degenerate;
    nit=RunData.IterCounter;
    npoll=RunData.PollSteps;
    spoll=RunData.SuccPollSteps;
    hits=RunData.Cache.hits;
    csize=RunData.Cache.size;
%catch
%end

display(x);
display(fx);