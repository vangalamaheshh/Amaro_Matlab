function c = dRanger_calculate_weirdpair_cutoff_from_isz(X,log10_cutoff)

require_fields(X,{'dat'});

if log10_cutoff>-2, error('log10_cutoff should be -2 or less'); end

s = sum(X.dat,2);
itcs = cumsum(s,1);
itcsf = bsxfun(@rdivide,itcs,itcs(end,:));
zz = log10(1-itcsf);
c = find(zz<=log10_cutoff,1);

if isinf(c), error('log10_cutoff is too stringent to be satisfied given isz distrib'); end

