function M = new_load_mutdata(maf,P)

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'medcovfile','/xchip/cga/reference/mutsig_params/mel_luad_lusc_mediancov.mat');
P = impose_default_value(P,'bagelfile','/xchip/cga/reference/mutsig_params/bagels.v2.mat');
P = impose_default_value(P,'mutation_blacklist','');
P = impose_default_value(P,'gene_min_frac_coverage_required',0.10);
P = impose_default_value(P,'consolidate_adjacent_muts',true);
P = impose_default_value(P,'consolidate_adjacent_muts_threshold',0); % only collapse mutations that are actually duplicates (dist=0)
P = impose_default_value(P,'apply_mutation_blacklist',true);

demand_file(maf);
demand_file(P.medcovfile);

fprintf('LOADING DATA...\n');

M=[];
if isfield(P,'isetname'), M.isetname = P.isetname; end
if isfield(P,'build'), M.build = P.build; end

% mutations
M.mut = load_struct(maf);
M.mut = add_and_convert_simple_fieldnames(M.mut);
%[M.mut P] = make_barebones(M.mut,P);
M.mut.gene = regexprep(M.mut.gene,'^([^\|]+).*','$1');
M.mut = add_helper_is_fields(M.mut);

if P.consolidate_adjacent_muts
  M.mut = collapse_adjacent_mutations(M.mut,P);
end

% patients
[M.patient.name tmp M.mut.pat_idx] = unique(M.mut.patient);
M.np = slength(M.patient);

% ttypes
if isfield(M.mut,'ttype')
  [M.ttype.name ui uj] = unique(M.mut.ttype);
  M.patient.ttype_idx = nansub(uj,listmap(M.patient.name,M.mut.patient));
  M.ttype.npat = histc(M.patient.ttype_idx,1:slength(M.ttype));
  M.ttype.nmut = histc(uj,1:slength(M.ttype));
  for i=1:slength(M.ttype)
    M.ttype.nmutgenes(i,1) = length(unique(M.mut.gene(nansub(M.patient.ttype_idx,M.mut.pat_idx)==i)));
  end
end

% mutation blacklist
if P.apply_mutation_blacklist
  try
    fprintf('Applying blacklist: %s\n',P.mutation_blacklist');
    M.mut = apply_mutation_blacklist(M.mut,P.mutation_blacklist,P);
  catch me
    fprintf('Error applying blacklist.\n');
    keyboard
  end
end

% impute callschemes
M = impute_callschemes(M);
count(M.patient.callscheme_name)
if mean(M.patient.callscheme>0)<0.1
  fprintf('WARNING:  No flanking mutation data available:  MutSigS2N/CV performance will suffer.\n');
end
if any(M.patient.callscheme==3)
  fprintf('WARNING:  Data appears to include WGS calls.  In future add WGS coverage model.\n');
end

% genes + territory + coverage
load(P.medcovfile,'D');
% remove genes with extremely low coverage (ideally should remove these from the medcovfile permanently)
frac = bsxfun(@rdivide,sum(D.gene_non_cov,3),sum(D.gene_non_terr,2));
frac = min(frac,[],2);
idx = find(frac<P.gene_min_frac_coverage_required);
if ~isempty(idx)
  oldng = D.ng;
  fprintf('Omitting %d/%d genes because they have extremely low coverage.\n',length(idx),oldng);
  D.gene = reorder_struct_exclude(D.gene,idx);
  D.ng = D.ng - length(idx);
  flds = fieldnames(D);
  for i=1:length(flds)
    sz = size(D.(flds{i}));
    if length(sz)>10, error('too many dimensions'); end
    if sz(1)==oldng, D.(flds{i})(idx,:,:,:,:,:,:,:,:,:)=[]; end
  end
end
M.cov = D; clear D;
M.gene = M.cov.gene;
M.ng = slength(M.gene);
M.mut.gene_idx = listmap(M.mut.gene,M.gene.name);

% bagels
if ~isempty(P.bagelfile) && ~strcmpi(P.bagelfile,'none')
  fprintf('Loading bagels.\n');
  bagels = load(P.bagelfile,'gene');
  map1 = listmap(M.gene.name,bagels.gene.name);
  map2 = listmap(bagels.gene.name,M.gene.name);
  M.gene.bagel = nansub(bagels.gene.bagel,map1);
  for i=1:size(M.gene.bagel,2), M.gene.bagel(:,i) = nansub(map2,M.gene.bagel(:,i)); end
  M.gene.nnei = sum(~isnan(M.gene.bagel),2);
end

% categories
if isfield(M.mut,'context65')
  if any(strcmp('',M.mut.context65))
    fprintf('Blanks in context65: will reload\n');
    M.mut = rmfield(M.mut,'context65');
  else
    M.mut = make_numeric(M.mut,'context65');
    bad = (M.mut.context65<1 | M.mut.context65>65 | isnan(M.mut.context65));
    if mean(bad)>0.1
      fprintf('Too many bad values in context65: will reload\n');
      M.mut = rmfield(M.mut,'context65');
    end
  end
end
if ~isfield(M.mut,'context65')
  if ~isfield(P,'build') || isempty(P.build)
    fprintf('Assuming hg19\n');
    P.build = 'hg19';
  end
  dirname = ['/cga/tcga-gsc/home/lawrence/db/' P.build '/context65'];
  fprintf('Loading context65 data from %s\n',dirname);
  M.mut.context65 = get_context(M.mut.chr,M.mut.start,dirname,P);
end
M.categ = M.cov.categ;
[M.mut.categ M.mut.categ_ignoring_null_categ] = assign_mut_categs(M.mut,M.categ);
M.mutclass = [M.categ.name' 'Total'];
M.ncat = slength(M.categ);
M.TOT = M.ncat+1;
M.NUM_INDEL_CLASSES = sum(grepmi('indel|null',M.categ.name));

% remove mutations outside gene/patient/categ set
bad = find(~(M.mut.categ>0 & M.mut.categ_ignoring_null_categ>0 & M.mut.pat_idx>0 & M.mut.gene_idx>0));
if ~isempty(bad)
  fprintf('Removing %d/%d mutations outside patient/gene/categ set\n',length(bad),slength(M.mut));
  M.mut = reorder_struct_exclude(M.mut,bad);
end

% mutation tallies and coverage tables
M = expand_to_n_and_N(M,P);




