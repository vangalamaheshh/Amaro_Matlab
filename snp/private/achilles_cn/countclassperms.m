function [count N Nt Nc] = countclassperms(classes,gtreat,gcontrol,binomial)
%COUNTCLASSPERMS count the number of permutations given class structure
%
% Given a partition of a set of scores into classes and selection of
% treatment and control groups, determine the theoretical number of distinct
% permutations of the scores.
%
%    [COUNT N Nt Nc] = countclassperms(CLASSES,GTREAT,GCONTROL,BINOMIAL)
%
% Returns the number of permutations in COUNT. The input argument CLASSES is a
% cell array containing vectors of score indices in each class. All indices
% in CLASSES should be distinct and their union should index all scores. 
% GTEST and GCONTROL are logical vectors defining treatment and control
% groups respectively. Their lengths equal the number of scores and their
% elements are true if the corresponding score belongs in the group. Scores
% should not belong to both groups, but some may belong to neither.
%
% The optional output N is a vector of counts of the number of scores in 
% each class and has the same length as CLASSES. Similarly, optional 
% outputs Nt and Nc count the number of treatment and control scores
% respectively.
%
if ~exist('binomial','var')
    binomial = true;
end
Nclasses = length(classes);
%% count permutations for each class
N = zeros(1,Nclasses);   % number of samples for each class 
Nt = zeros(1,Nclasses);  % number of treatment samples for each class
Nc = zeros(1,Nclasses);  % number of control samples for each class
% calculate permutation count to decide exact vs Monte Carlo
count = 1;
for c = 1:Nclasses
    % each class multiplies in a trinomial coefficient's worth
    % of permutations
    N(c) = length(classes{c});
    Nt(c) = sum(gtreat(classes{c}));
    Nc(c) = sum(gcontrol(classes{c}));
    if binomial
        clscount = factorial(Nt(c)+Nc(c)) / (factorial(Nt(c))*factorial(Nc(c)));
    else % trinomial
        clscount = factorial(N(c)) / ...
            ( factorial(Nt(c)) * factorial(Nc(c)) * factorial(N(c)-Nt(c)-Nc(c)) );
    end
    count = count * clscount;
end
