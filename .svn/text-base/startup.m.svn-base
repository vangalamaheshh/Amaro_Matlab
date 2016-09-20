function startup
% STARTUP adds CancerGenomeAnalysis directories to path on matlab startup.

%
if ~isdeployed
    disp([mfilename('fullpath') ': Adding CancerGenomeAnalysis directories to matlab path']);

    p = genpath('~/CancerGenomeAnalysis/trunk/matlab/');

    p = textscan(p,'%[^:]','Delimiter',':');
    bool = cellfun(@isempty,(regexp(p{:},'.svn')));
    p = p{1}(bool);
    p = strcat(p,':');
    addpath(cell2mat(p'));
end


