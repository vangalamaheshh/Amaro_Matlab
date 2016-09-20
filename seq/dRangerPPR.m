function dRangerPPR(sample,P)
if ~exist('P','var'),P=[]; end
dRangerPreprocess(sample,P);
dRangerPrepare(sample,P);
dRangerRun(sample,P);
