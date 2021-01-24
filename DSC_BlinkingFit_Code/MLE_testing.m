% MLE_testing.m
%
% Script for ascertaining how much of a difference MLE might actually make
% in our fits
%
% Nathanael Kazmierczak, 04/10/2020

s1 = [0,0,5];
s2 = [1,2,2];

loglik = @(lam,ks) sum(ks)*log(lam) - 3*lam - sum(factorial(ks));
loglik1 = @(lam) -1*loglik(lam,s1);
loglik2 = @(lam) -1*loglik(lam,s2);

best_lambda1 = fminsearch(loglik1,2)
best_lambda2 = fminsearch(loglik2,1)
