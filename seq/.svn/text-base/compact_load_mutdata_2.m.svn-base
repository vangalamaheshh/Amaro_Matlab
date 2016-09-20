function M = compact_load_mutdata_2(maf,P)
% replacement for new_load_mutdata_2
%
% --> no more large n and N tables!
% --> represents N as an indexed form in .cov, with pointers from .pat
% --> keeps all mutation data in .mut, does not generate a full n table.
%
% using new system of coverage (ncd/syn/mis/non/spl)

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'covmodelfile','/xchip/cga/reference/mutsig_params/meancov.v4.mat');
P = impose_default_value(P,'catfile',['/cga/tcga-gsc/home/lawrence/cga/trunk/matlab/seq/'...
                    'categs_CpGtransit_CpGtransver_otherCGtransit_otherCGtransver_ATtransit_ATtransver.txt']);
P = impose_default_value(P,'covarfile','/cga/tcga-gsc/home/lawrence/hn/analysis/20121004/covars.with_transform.v3.mat');
P = impose_default_value(P,'covars_to_use',22:26); % transformed: expr, RT, gc, hiC, gdens1mb
P = impose_default_value(P,'gene_min_frac_coverage_required',0.10);
P = impose_default_value(P,'consolidate_adjacent_muts',true);
P = impose_default_value(P,'consolidate_adjacent_muts_threshold',0); % only collapse mutations that are actually duplicates (dist=0)
P = impose_default_value(P,'apply_mutation_blacklist',true);
P = impose_default_value(P,'mutation_blacklist','');

demand_file(maf);
demand_file(P.covmodelfile);
demand_file(P.covarfile);
demand_file(P.catfile);

fprintf('LOADING DATA\n');

M=[];
if isfield(P,'isetname'), M.isetname = P.isetname; end
if isfield(P,'build'), M.build = P.build; end

% mutations
fprintf('Loading mutations...\n');
M.mut = load_struct(maf);
fprintf('Converting mutation data...\n');
M.mut = add_and_convert_simple_fieldnames(M.mut);
M.mut.gene = regexprep(M.mut.gene,'^([^\|]+).*','$1'); % legacy fix
fprintf('Adding helper fields...\n');
M.mut = add_helper_is_fields(M.mut);
if P.consolidate_adjacent_muts, M.mut = collapse_adjacent_mutations(M.mut,P); end

% calculate M.mut.effect
M.mut.effect = repmat({'---'},slength(M.mut),1);
idx = find(M.mut.is_flank); M.mut.effect(idx) = repmat({'ncd'},length(idx),1);
idx = find(M.mut.is_silent); M.mut.effect(idx) = repmat({'syn'},length(idx),1);
idx = find(M.mut.is_missense); M.mut.effect(idx) = repmat({'mis'},length(idx),1);
idx = find(M.mut.is_nonsense); M.mut.effect(idx) = repmat({'non'},length(idx),1);
idx = find(M.mut.is_splice); M.mut.effect(idx) = repmat({'spl'},length(idx),1);
% the following miscellaneous other mutation types are also treated like splice-site mutations:
miscnull = grepi('non.?stop|read.?thr|start',M.mut.type,1);
M.mut.effect(miscnull) = repmat({'spl'},length(miscnull),1);
% indel types
idx = find(M.mut.is_indel & ~M.mut.is_coding); M.mut.effect(idx) = repmat({'indel_ncd'},length(idx),1);
idx = find(M.mut.is_indel & M.mut.is_coding); M.mut.effect(idx) = repmat({'indel_cod'},length(idx),1);
idx = find(M.mut.is_indel & M.mut.is_splice); M.mut.effect(idx) = repmat({'indel_spl'},length(idx),1);
allowed_effects = {'ncd','syn','mis','non','spl','indel_ncd','indel_spl','indel_cod'};
idx = find(~ismember(M.mut.effect,allowed_effects));
if ~isempty(idx)
  fprintf('Note: removing the following unhandled-type mutations:\n');
  count(M.mut.type(idx),1)
  M.mut = reorder_struct_exclude(M.mut,idx);
end

% patients
[M.pat.name tmp M.mut.pat_idx] = unique(M.mut.patient);
M.np = slength(M.pat);

% ttypes
if isfield(M.mut,'ttype')
  [M.ttype.name ui uj] = unique(M.mut.ttype);
  M.pat.ttype_idx = nansub(uj,listmap(M.pat.name,M.mut.patient));
  M.pat.ttype = nansub(M.ttype.name,M.pat.ttype_idx);
  M.ttype.npat = histc(M.pat.ttype_idx,1:slength(M.ttype));
  M.ttype.nmut = histc(uj,1:slength(M.ttype));
  for i=1:slength(M.ttype)
    M.ttype.nmutgenes(i,1) = length(unique(M.mut.gene(nansub(M.pat.ttype_idx,M.mut.pat_idx)==i)));
  end
  M.nttype = slength(M.ttype);
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
count(M.pat.callscheme_name)
if mean(M.pat.callscheme>0)<0.1
  fprintf('WARNING:  No flanking mutation data available:  MutSigS2N/CV performance will suffer.\n');
end

% genes + coverage
fprintf('Loading coverage models...\n');
load(P.covmodelfile,'C');
demand_fields(C,{'type','ntype','gene','cat','ncat','gene_effect_terr','gene_effect_cov'});
if C.ncat~=1885, error('wrong C.ncat'); end
if any(isinf(C.gene_effect_terr(:))) || any(isinf(C.gene_effect_cov(:))), error('infs in terr/cov'); end
if any(isnan(C.gene_effect_terr(:))) || any(isnan(C.gene_effect_cov(:))), error('nans in terr/cov'); end
M.cov = C; clear C;

% coverage indexing to models
if M.cov.ntype~=5, error('wrong number of coverage models'); end
M.pat.cov_idx = nan(M.np,1);
M.pat.cov_idx(M.pat.callscheme==0) = 1;         % coding only
M.pat.cov_idx(M.pat.callscheme==1) = 2;         % exome with typical interval list (~100bp flanks)
M.pat.cov_idx(M.pat.callscheme==2) = 3;         % all capture (no interval list)
M.pat.cov_idx(M.pat.callscheme==3) = 5;         % all genome (WGS)
% TO DO: add exomeplus as a detected callscheme

% remove genes with very low coverage of coding regions
% (based on the coverage models that were identified as being relevant)
cod = grepv('noncoding',M.cov.cat.name,1);
gene_terr = sum(M.cov.gene_effect_terr(:,:,cod),3);
gene_model_cov = sum(M.cov.gene_effect_cov(:,:,cod),3);
num_each_model = histc(M.pat.cov_idx,(1:M.cov.ntype)');
gene_tot_terr = gene_terr * M.np;
gene_tot_cov = gene_model_cov * num_each_model;
gene_frac_cov = gene_tot_cov ./ gene_tot_terr;
idx = find(gene_frac_cov<P.gene_min_frac_coverage_required);
if ~isempty(idx)
  oldng = M.cov.ng;
  fprintf('Omitting %d/%d genes because they have extremely low coverage.\n',length(idx),oldng);
  M.cov.gene = reorder_struct_exclude(M.cov.gene,idx);
  M.cov.ng = M.cov.ng - length(idx);
  M.cov.gene_effect_cov(idx,:,:)=[];
  M.cov.gene_effect_terr(idx,:,:)=[];
end

% final genelist
M.gene = M.cov.gene; M.cov = rmfield(M.cov,'gene');
M.ng = slength(M.gene); M.cov = rmfield_if_exist(M.cov,'ng');
M.mut.gene_idx = listmap(M.mut.gene,M.gene.name);

% categories (category list is now provided as an input; can collapse coverage categories however we like)
M.categ = load_struct(P.catfile);
M.ncat = slength(M.categ);
if ~isempty(grepi('null|indel',M.categ.name)), error('catfile should not contain indel/null category'); end
if isfield(M.categ,'type') && ~all(strcmp('point',M.categ.type)), error('catfile should contain only point categories'); end

% collapse coverage model to (categ x 5degree) table
% two steps: 1. collapse 1885 to 192x5degree
%            2. collapse 192x5degree to ncatx5degree
M.cov.gene_effect_terr = collapse_1885_to_192x5(M.cov.gene_effect_terr,3);
M.cov.gene_effect_terr = collapse_192_to_categ_set(M.cov.gene_effect_terr,M.categ,3);
M.cov.gene_effect_cov = collapse_1885_to_192x5(M.cov.gene_effect_cov,3);
M.cov.gene_effect_cov = collapse_192_to_categ_set(M.cov.gene_effect_cov,M.categ,3);
M.cov.dims_of_gene_effect_cov = 'gene, covmodel, ncat, effect (ncd/syn/mis/non/spl)';
M.cov = move_field_to_before(M.cov,'dims_of_gene_effect_cov','gene_effect_cov');
M.cov = rmfield(M.cov,{'cat','ncat'}); % remove old category information

% assign mutations to categories
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
end,end,end
if ~isfield(M.mut,'context65')
  if ~isfield(P,'build') || isempty(P.build)
    fprintf('Assuming hg19\n');
    P.build = 'hg19';
  end
  dirname = ['/cga/tcga-gsc/home/lawrence/db/' P.build '/context65'];
  fprintf('Loading context65 data from %s\n',dirname);
  M.mut.context65 = get_context(M.mut.chr,M.mut.start,dirname,P);
end
M.mut = rmfield_if_exist(M.mut,{'categ','categ_ignoring_null_categ'});
fprintf('Note: Indels will be handled separately from "categ" system (ignore next warning about indels).\n');
[tmp M.mut.categ] = assign_mut_categs(M.mut,M.categ);
% (everything downstream ignores any "null category")

% remove mutations outside gene/patient/categ set
bad = find(~((M.mut.categ>0 | M.mut.is_indel) & M.mut.pat_idx>0 & M.mut.gene_idx>0));
if ~isempty(bad)
  fprintf('Removing %d/%d mutations outside patient/gene/categ set\n',length(bad),slength(M.mut));
  M.mut = reorder_struct_exclude(M.mut,bad);
end

% covariates
load(P.covarfile,'V');
idx = listmap(M.gene.name,V.gene_names); V=rmfield(V,'gene_names');
for i=1:length(V.val), V.val{i} = nansub(V.val{i},idx); end; M.V = reorder_struct(V,P.covars_to_use);

%%%
flds = {'build','np','nttype','ng','ncat','mut','pat','ttype','cov','gene','categ','V'};
flds(~ismember(flds,fieldnames(M)))=[];
M = order_fields(M,flds);






