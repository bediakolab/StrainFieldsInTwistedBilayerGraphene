% testTrigFittingFunction.m
%
% Nathanael Kazmierczak, 12/28/2019

testfun = @(v)blinkingFittingFunction(v) - trigFittingFunctions(v);
testfun([0.9,-0.38])
testfun([-1,-0.12])
testfun([1.3,0])
% looks good 
