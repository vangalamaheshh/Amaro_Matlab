function N = analyze_territory(M,P)
%
% given M.targ, with chr,strand,start,end,frame
% computes matrix N (ntarg,12,16,2)
% where dim1 is target no
%       dim2 is mutation types (e.g. A->C)
%       dim3 is mutation contexts (e.g. A_A)
%       dim4 is 1=silent 2=nonsilent
%
% Mike Lawrence 2008-09-08
%
if ~exist('P','var'), P=[]; end
P=impose_default_value(P,'targ_genome_build','hg18');

build = P.targ_genome_build;
require_fields(M,'targ');
target = M.targ;
ntarg = length(target.gene);

N = zeros(ntarg,12,16,2);    % target x mutation-type x base-context x (silent/nonsilent)
bases='ACGT';
tic
for t=1:ntarg
  if ~mod(t,100), a=toc; fprintf('%d/%d ', t,ntarg); end

  seq = genome_region(target.chr(t), target.start(t)-3, target.end(t)+3, build);

  if strcmp(target.strand{t}, '-'), seq = my_seqrcomplement(seq); end

  % find DNA contexts
  context = survey_contexts(seq,1);   % context_type=1

  % find change consequences
  codon_start = target.frame(t);
  for i=3:length(seq)-3
    if i==codon_start+3, codon_start = i; end
    if i<codon_start || codon_start+2 > length(seq), continue; end
    x = find(seq(i)==bases);
    col = context(i);
    newbases = setdiff(bases,seq(i));
    old_codon = seq(codon_start:codon_start+2);
    old_aa = my_nt2aa(old_codon);
    for j=1:3
      newbase = newbases(j);
      row = (x-1)*3+j;
      new_dna = seq;
      new_dna(i) = newbase;
      new_codon = new_dna(codon_start:codon_start+2);
      new_aa = my_nt2aa(new_codon);
      silent = strcmp(old_aa,new_aa);
      if silent, page=1; else page=2; end
      N(t,row,col,page)=N(t,row,col,page)+1;
    end % next j
  end % next i
end   % next target
