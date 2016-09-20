function m = fix_splice(m)
% any mutations within +-2bp (-1 0 +1 +2) of a splice site must be called splice_site!

if ~isfield(m,'splicedist'), m.splicedist = get_splicedist(m.chr,m.pos); end

midx=find(m.splicedist>=-1 & m.splicedist<=2 & ~grepmi('splice',m.type));
fprintf('Fixing %d splice-site mutations\n',length(midx));
m.type(midx) = repmat({'Splice_site'},length(midx),1);

