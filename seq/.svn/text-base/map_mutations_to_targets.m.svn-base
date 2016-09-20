function map = map_mutations_to_targets(x,T,P)

if ~exist('P','var'), P=[]; end

fprintf('Mapping mutations to targets: ');

T.chr = convert_chr(T.chr,P);
T = make_numeric(T,{'start','end'});
x.chr = convert_chr(x.chr,P);
x = make_numeric(x,{'start','end'});

map = nan(slength(x),1);
for chr=1:get_chrct(P),fprintf('chr%d ',chr);
  tidx = find(T.chr==chr);
  xidx = find(x.chr==chr);
  for i=1:length(xidx), j = xidx(i);
    idx = tidx(find(x.start(j)<=T.end(tidx) & x.end(j)>=T.start(tidx),1));
    if ~isempty(idx), map(j) = idx; end
  end
end, fprintf('\n');
