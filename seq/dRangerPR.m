function dRangerPR(sample,P)
if ~exist('P','var'),P=[]; end
dRangerPreprocess(sample,P);
dRangerRun(sample,P);
