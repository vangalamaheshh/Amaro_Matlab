function M = impute_frames(M,P)

if ~exist('P','var'), P=[]; end
P=impose_default_value(P,'targ_genome_build','hg18');

build = P.targ_genome_build;
require_fields(M,'targ');
target = M.targ;
ntarg = length(target.gene);
target.genename = M.gene.name(target.gene);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   BEGIN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Imputing frames...\n');
frame = zeros(ntarg,1);

% first pass: compare to proteome fasta file

fprintf('Loading proteome fasta file\n');
P = load_fasta('/xchip/tcga/gbm/analysis/lawrence/db/human_proteins.fa');
P.gene = regexprep(P.header,'[\(\)\>]','');

pg = unique(P.gene);
for t=1:ntarg
  if isnan(target.chr(t)), continue; end
  if ~mod(t,10), fprintf('%d/%d ', t,ntarg); end
  result = 0;
  pno = grep(['^' target.genename{t} '$'],P.gene,1);
  for p=1:length(pno)
    seq = genome_region(target.chr(t), target.start(t), target.end(t), build);
    if strcmp(target.strand{t},'-'), seq = my_seqrcomplement(seq); end
    for f=1:3
      aa{f} = my_nt2aa(seq(f:end));
      if ~isempty(aa{f})
        tmp = grep(aa{f},P.seq{pno(p)},1);
        if ~isempty(tmp), result=f; break; end
      end
    end
  end
  frame(t) = result;
end
fprintf('\n');


M.targ.frame = frame;

