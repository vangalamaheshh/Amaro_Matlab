function find_clustered_mutations(M,P,R)
% Mike Lawrence 2010

if ~exist('P','var'),P=[]; end
P = impose_default_value(P,'build','hg18');
P = impose_default_value(P,'build_refseq',P.build);
P = impose_default_value(P,'build_genome_region',P.build);
P = impose_default_value(P,'clustered_mutations_output','*required*');
P = impose_default_value(P,'method',1);

if ~exist('R','var'), R = load_refseq(P.build_refseq); end
require_field(R,'tx_len');

if P.method==1
  if isfield(M,'use_nonsilent')
    mut = reorder_struct(M.mut,M.use_nonsilent);
  else
    mut = M.mut;
    mut = reorder_struct(mut,grepvi('Sil|Syn',mut.type,1));
  end
  mut = make_numeric(mut,{'start','end'});
  
  if ~isfield(mut,'change')
    if isfield(mut,'tum_allele1') && isfield(mut,'tum_allele2')
      mut.change = find_change_from_ref_tum1_tum2(mut.start,mut.end,mut.ref_allele,mut.tum_allele1,mut.tum_allele2);
    else
      mut.change = find_change_from_ref_tum1_tum2(mut.start,mut.end,mut.ref_allele,mut.newbase,mut.newbase);
    end
    idx = grepi('error',mut.change,1);
    if ~isempty(idx)
      fprintf('Ignoring the mutations listed above.\n');
      mut = reorder_struct_exclude(mut,idx);
    end
  end
  
  if isnumeric(mut.gene)
    mut.gene = nansub(M.gene.name,mut.gene);
  end
  
  [g gi gj] = unique(mut.gene);
  
  X = [];
  X.num = (1:length(g))';
  X.gene = g;
  X.desc = get_longnames(g);
  X.n = nan(length(g),1);
  X.mindist = nan(length(g),1);
  X.nmuts0 = nan(length(g),1);
  X.nmuts3 = nan(length(g),1);
  X.nmuts12 = nan(length(g),1);
  X.npairs0 = nan(length(g),1);
  X.npairs3 = nan(length(g),1);
  X.npairs12 = nan(length(g),1);
  
  fprintf('Analyzing clustered mutations\n');
  PP = [];
  PP.build = P.build;
  PP.build_genome_region = P.build_genome_region;
  PP.include_new_fields = true;
  for i=1:length(g), if ~mod(i,100), fprintf('%d/%d ',i,length(g)); end
    midx = find(gj==i);
    if length(midx)<2, continue; end
    muti = reorder_struct(mut,midx);
    chose_transcript = false;
    if isfield(muti,'transcript') % if transcript annotations available, choose most-used transcript
      [t ti tj] = unique(muti.transcript);
      [tmp0,w] = max(histc(tj,1:length(t)));
      ridx = grep(t{w},R.transcript,1);
      if length(ridx)>=1
        ridx=ridx(1);
        chose_transcript = true;
      end
    end
    if ~chose_transcript   % else choose refseq transcript with longest sum of exon lengths
      rgidx = find(strcmpi(g{i},R.gene));
      if ~isempty(rgidx)
        [tmp,idx] = max(R.tx_len(rgidx));
        ridx = rgidx(idx);
        chose_transcript = true;
      end
    end
    if ~chose_transcript, continue; end
    Rw = reorder_struct(R,ridx);
    len = compute_protein_length(Rw.code_start,Rw.code_end,Rw.exon_starts{1},Rw.exon_ends{1});
    if slength(muti)>=100, fprintf('<%s muts: ',g{i}); end
    try
      muti = classify_muts(muti,PP,Rw);
    catch me
      fprintf('Problem with re-annotating mutations in gene %s: skipping.\n',g{i});
      continue;
    end
    if slength(muti)>=100, fprintf('/> '); end
    d = dist(muti.proteinpos,[]);
    for k=1:size(d,1), d(k,k)=inf; end
    X.n(i) = length(midx);
    X.mindist(i) = min(d(:));
    X.nmuts0(i) = sum(sum(d==0,2))/2;
    X.nmuts3(i) = sum(sum(d<=3,2))/2;
    X.nmuts12(i) = sum(sum(d<=12,2))/2;
    X.npairs0(i) = sum(d(:)==0)/2;
    X.npairs3(i) = sum(d(:)<=3)/2;
    X.npairs12(i) = sum(d(:)<=12)/2;
  end, fprintf('\n');
  X = sort_struct(X,{'mindist','nmuts0','nmuts3','nmuts12','npairs0','npairs3','npairs12'},[1 -1 -1 -1 -1 -1 -1]);
  X = reorder_struct(X,X.mindist<inf);
  save_struct(X,P.clustered_mutations_output);

elseif P.method==2

  % new method
  error('new method not written yet');

else
  error('unknown P.method');
end


