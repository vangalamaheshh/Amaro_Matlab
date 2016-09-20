function a = fixlog10(a,infvalue)

if ~exist('infvalue','var'), infvalue = -1; end

idx0 = find(a==0);

a = log10(a);

a(idx0) = infvalue;

