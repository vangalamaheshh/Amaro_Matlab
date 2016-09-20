function cDNA_pos = find_cDNA_pos(M,T)

require_fields(M,{'chr','start'});
M.pos = M.start;
require_fields(T,{'chr','start','end'});

if ~issorted(T.start), error('T needs to be sorted'); end
nm = slength(M);
nt = slength(T);

M.chr = convert_chr(M.chr);
M = make_numeric(M,{'start','end'});
T.chr = convert_chr(T.chr);
T = make_numeric(T,{'start','end'});

T.len = T.end-T.start+1;
T.clen = cumsum(T.len);

cDNA_pos = nan(nm,1);
for i=1:nm
  tidx = find(M.pos(i)>=T.start & M.pos(i)<=T.end & M.chr(i)==T.chr,1);
  cDNA_pos(i) = T.clen(tidx) - (T.end(tidx) - M.pos(i));

%  fprintf('Mutation at %d is in exon %d (%d - %d) [%d - %d], cDNA_pos = %d\n',...
%    M.pos(i),tidx,T.start(tidx),T.end(tidx),T.clen(tidx)-T.len(tidx)+1,T.clen(tidx),cDNA_pos(i));
end

