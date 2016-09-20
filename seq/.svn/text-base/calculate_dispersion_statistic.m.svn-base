function calculate_dispersion_statistic(M,P)
% calculate_dispersion_statistic(M,P)
%
% Procedure:
%   1. consider all nonsynonymous coding mutations   
%   2. consider all genes with more than one patient having such a mutation
%
%
% input data used:
%     M.mut  (patient,type,transcript,proteinchange)
%
%
% Mike Lawrence 2010
if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'same_sample_single_event_threshold',6);
P = impose_default_value(P,'build','hg18');

R = load_refseq(P.build);

% m = M.mut;

m = load_struct(['/xchip/cga1/firehose_output/Individual_Set/An_OV_Capture_TCGA_123/capture/mut/' ...
                 'An_OV_Capture_TCGA_123.maf.annotated']);

m.chr = convert_chr(m.chr);
m = reorder_struct(m,~isnan(m.chr));
m.pos = str2double(m.start);
m.newbase = m.tum_allele1;
idx = find(strcmp(m.ref_allele,m.tum_allele1));
m.newbase(idx) = m.tum_allele2(idx);

[g gi gj] = unique(m.gene);
for i=1:length(g)           %     i=grep('^TP53$',g,1)
  fprintf('%s\n',g{i});
  midx = find(gj==i);
  if length(midx)<2, continue; end
  % pick majority-vote transcript
  [t ti tj] = unique(m.transcript(midx));
  [tmp0,w] = max(histc(tj,1:length(t)));
  ridx = grep(t{w},R.transcript,1);
  if isempty(ridx), error('%s no longer found in Refseq',t{w}); end
  if length(ridx)>1, error('%s has multiple entries in Refseq',t{w}); end
  Rw = reorder_struct(R,ridx);
  len = compute_protein_length(Rw.code_start,Rw.code_end,Rw.exon_starts{1},Rw.exon_ends{1});
  % reannotate mutations if necessary
  if length(t)>1
    nw = midx(tj~=w);
    mi = reorder_struct(m,nw);
    mi = classify_muts(mi,struct('include_new_fields',true),Rw);
    m = struct_assign(m,nw,mi);
  end
  % keep only missense, nonsense, splice-site
  

  % collapse nearby same-sample events
  

end



