clear Problem;
clear Options;
clear InitPop;

% Problem definition
Problem.ObjFunction='hs024_obj';
Problem.A = [-1/sqrt(3) 1;
    -1 -sqrt(3);
    1 sqrt(3)];
Problem.b = [0;0;6];

Problem.LB = [0;0];
Problem.UB = [Inf;Inf];

% Initial guess

InitPop(1).x=[3; 0.5];
InitPop(2).x=[1; 0.5];

% Algorithm options
Options.Size=40;
%Options.Size=1;
Options.MaxIter=2000;
Options.MaxObj=2000;
%Options.SearchType=0;
Options.LoadCache=0;

% Run the algorithm
[x,fx,RunData]=PSwarm(Problem, InitPop, Options);
nfo=RunData.ObjFunCounter
deg=RunData.Degenerate
x
fx