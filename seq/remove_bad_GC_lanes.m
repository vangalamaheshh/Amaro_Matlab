function B = remove_bad_GC_lanes(B)
% Mike Lawrence
% 2009-07-30

B.raw = B.raw(:,B.good,:,:);
B.dat = B.dat(:,B.good,:);
B.cov = B.cov(B.good,:);
B.datnorm = B.datnorm(:,B.good,:);
B.lanes = reorder_struct(B.lanes,B.good);
B.good = (1:B.nlanes)';
B.nlanes = length(B.good);
