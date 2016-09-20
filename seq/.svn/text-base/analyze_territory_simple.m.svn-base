function N = analyze_territory_simple(M,P)
%
% given M.targ, with chr,start,end
% computes matrix N (ntarg,3)
% where dim1 is target no
%       dim2 is mutation type: (1) CpG (2) C+G (3) A+T
%
% Mike Lawrence 2008-09-30
%
if ~exist('P','var'), P=[]; end
P=impose_default_value(P,'targ_genome_build','hg18');

build = P.targ_genome_build;
require_fields(M,'targ');
target = M.targ;
ntarg = slength(target);

N = zeros(ntarg,3);
for t=1:ntarg
  if ~mod(t,1000), fprintf('%d/%d ', t,ntarg); end
  if isnan(target.chr(t)), continue; end
  seq = upper(genome_region(target.chr(t), target.start(t)-1, target.end(t)+1, build));
  seq = regexprep(seq,'CG','XX');
  seq = seq(2:end-1);
  N(t,1) = sum(seq=='X');
  N(t,2) = sum(seq=='C')+sum(seq=='G');
  N(t,3) = length(seq)-N(t,1)-N(t,2);
end


