function B = combine_GC_TN(T,N)
% Mike Lawrence
% 2009-07-30

B = T;
B.raw = cat(2,T.raw,N.raw);
B.dat = cat(2,T.dat,N.dat);
B.nlanes = T.nlanes + N.nlanes;
B.cov = [T.cov;N.cov];
B.good = [T.good; T.nlanes+N.good];
B.datnorm = cat(2,T.datnorm,N.datnorm);
T.lanes.tn = repmat({'t'},T.nlanes,1);
N.lanes.tn = repmat({'n'},N.nlanes,1);
if isfield(N.lanes,'col1'), N.lanes.col1 = N.lanes.col1 + T.nlanes; end
if isfield(N.lanes,'lane'), N.lanes.lane = N.lanes.lane + T.nlanes; end
B.lanes = combine_structs({T.lanes,N.lanes});
