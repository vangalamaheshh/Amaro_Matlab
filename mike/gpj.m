function f = gpj(u)
% f = gpj(u)
%
% example:
% gpj('http://cga-gdac.broadinstitute.org:8090/gp/jobResults/122192/PR-1783-Tumor.breakpoints.txt')
%
% ans = 
%      '/xchip/cga/gdac-prod/genepattern/jobResults/122192/PR-1783-Tumor.breakpoints.txt'
%
f = regexprep(u,'.*(jobResults/.*)$','/xchip/cga/gdac-prod/genepattern/$1');

if iscell(f)
  for i=1:length(f)
    if ~exist(f{i},'file'), fprintf('WARNING file not found: %s\n',f{i}); end
  end
else
  if ~exist(f,'file'), fprintf('WARNING file not found: %s\n',f); end
end
