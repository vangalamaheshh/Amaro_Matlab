function G = MutSigCV_standalone_code_6_1D(M,gta,P,outdir)
% uses COMPACT representations
% with 1D projection
% restores the use of bagels
% uses (ncd,syn,mis,non,spl) system

fprintf('%s [top]\n',datestr(now));

flaghash = sparse([]);

% omitted gta?
if nargin==2 && (isempty(gta)||isstruct(gta)), P=gta; clear gta; end
if nargin==3 && (isempty(gta)||isstruct(gta)), outdir=P; P=gta; clear gta; end

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'debug',false); DEBUG=P.debug;
P = impose_default_value(P,'collapse_categories',false);
if P.collapse_categories, error('not implemented'); end
if isfield(P,'categories_to_use'), error('not implemented'); end
P = impose_default_value(P,'ttypes_to_use',[]);
P = impose_default_value(P,'ttype_grep','');
P = impose_default_value(P,'ttype_grepv','');
P = impose_default_value(P,'patients_to_use',[]);
P = impose_default_value(P,'ignore_indels',false);
P = impose_default_value(P,'scaling',true);
P = impose_default_value(P,'impute_full_cov_when_promotes_significance',false);
P = impose_default_value(P,'increase_terr_to_max_cov_represented',true);
P = impose_default_value(P,'special_treatment_for_category',[]);
if ~isempty(P.special_treatment_for_category), error('not implemented'); end
P = impose_default_value(P,'use_bagels',true);
if ~P.use_bagels, error('not implemented'); end
P = impose_default_value(P,'use_bagels_from',[]);
P = impose_default_value(P,'per_category_bagels',false);
if P.per_category_bagels, error('not implemented'); end
P = impose_default_value(P,'remove_superfluous_signal',false);
P = impose_default_value(P,'min_neighbors',0);
P = impose_default_value(P,'max_neighbors',50);
P = impose_default_value(P,'min_territory_ratio',0);
P = impose_default_value(P,'qual_min',0.05);
P = impose_default_value(P,'indel_min_neighbors',1);
if P.indel_min_neighbors<1, error('indel_min_neighbors must be at least 1'); end
P = impose_default_value(P,'indel_min_territory_ratio',0);
P = impose_default_value(P,'indel_max_neighbors',50);
P = impose_default_value(P,'indel_qual_min',0.05);
P = impose_default_value(P,'num_neighbor_patients',0);
P = impose_default_value(P,'null_score_boost',0);
P = impose_default_value(P,'numbins_scheme',2);
P = impose_default_value(P,'signal','nonsilent');
P = impose_default_value(P,'score_multiplier',10);
if isfield(P,'minimum_effect_size') && ~isfield(P,'min_effect_size'), P=rename_field(P,'minimum_effect_size','min_effect_size'); end
P = impose_default_value(P,'min_effect_size',1);
if any(P.min_effect_size<1), error('P.min_effect_size must be 1.00 or greater'); end
P = impose_default_value(P,'forced_bagels',{});
P = impose_default_value(P,'output_P0_table',false);

if exist('outdir','var')
  ede(outdir);
  progress_f = fopen([outdir '/progress.txt'],'wt');
  notes_f = fopen([outdir '/notes.txt'],'wt');
end

if ~isfield(M,'pat') & isfield(M,'patient'), M=rename_field(M,'patient','pat'); end
if ~isfield(M,'ng'), M.ng = slength(M.gene); end
if ~isfield(M,'np'), M.np = slength(M.pat); end
if ~isfield(M,'ncat'), M.ncat = slength(M.categ); end
demand_fields(M,{'mut','pat','gene','ng','np','ncat','categ'});
demand_fields(M.cov,{'gene_effect_cov','ntype'});
demand_fields(M.pat,{'cov_idx'});
if ~isfield(M.mut,'effect')
  % calculate M.mut.effect
  M.mut.effect = repmat({'---'},slength(M.mut),1);
  idx = find(M.mut.is_flank); M.mut.effect(idx) = repmat({'ncd'},length(idx),1);
  idx = find(M.mut.is_silent); M.mut.effect(idx) = repmat({'syn'},length(idx),1);
  idx = find(M.mut.is_missense); M.mut.effect(idx) = repmat({'mis'},length(idx),1);
  idx = find(M.mut.is_nonsense); M.mut.effect(idx) = repmat({'non'},length(idx),1);
  idx = find(M.mut.is_splice); M.mut.effect(idx) = repmat({'spl'},length(idx),1);
  % the following miscellaneous other mutation types are also treated like splice-site mutations:
  miscnull = grepi('non.?stop|read.?thr|start',M.mut.type,1);   % NOTE: remove this "start" category for the next release
  M.mut.effect(miscnull) = repmat({'spl'},length(miscnull),1);
  % indel types
  idx = find(M.mut.is_indel & ~M.mut.is_coding); M.mut.effect(idx) = repmat({'indel_ncd'},length(idx),1);
  idx = find(M.mut.is_indel & M.mut.is_coding); M.mut.effect(idx) = repmat({'indel_cod'},length(idx),1);
  idx = find(M.mut.is_indel & M.mut.is_splice); M.mut.effect(idx) = repmat({'indel_spl'},length(idx),1);
end
demand_fields(M.mut,{'gene_idx','pat_idx','effect','is_indel'});
demand_fields(M.mut,{'chr','start'});  % for finding nsite
allowed_effects = {'ncd','syn','mis','non','spl','indel_ncd','indel_spl','indel_cod'};
idx = find(~ismember(M.mut.effect,allowed_effects));
if ~isempty(idx)
  fprintf('Note: removing the following unhandled-type mutations:\n');
  count(M.mut.type(idx),1)
  M.mut = reorder_struct_exclude(M.mut,idx);
end
if isfield(M.mut,'categ_ignoring_null_categ')
  M.mut.categ = M.mut.categ_ignoring_null_categ;
end

sz = size(M.cov.gene_effect_cov);
if length(sz)~=4 || ~all(sz==[M.ng M.cov.ntype M.ncat 5]), error('M.cov.gene_effect_cov is wrong size'); end

fprintf('%s [select]\n',datestr(now));

if ~isempty(P.ttypes_to_use) && ~isempty(P.patients_to_use), error('please use only one'); end

% which ttypes to use
if ~isempty(P.ttypes_to_use) && (~isempty(P.ttype_grep)||~isempty(P.ttype_grepv)), error('please use only one'); end
if ~isempty(P.ttype_grep)||~isempty(P.ttype_grepv)
  orig_npat = slength(M.pat); pidx = (1:orig_npat)'; M.pat.idx_orig = pidx;
  if ~isempty(P.ttype_grep), pidx = pidx(grepm(P.ttype_grep,M.pat.ttype(pidx))); end
  if ~isempty(P.ttype_grepv), pidx = pidx(~grepm(P.ttype_grepv,M.mut.ttype(pidx))); end
  M.pat = reorder_struct(M.pat,pidx); M.np = length(pidx);
  map = nan(orig_npat,1); map(M.pat.idx_orig) = (1:slength(M.pat));
  M.mut.pat_idx = map(M.mut.pat_idx);
  M.mut = reorder_struct_exclude(M.mut,isnan(M.mut.pat_idx));
  notes_fprintf('Using specified grep/grepv on ttypes (%d patients).\n',length(pidx));
elseif ~isempty(P.ttypes_to_use)
  if ~isnumeric(P.ttypes_to_use), error('P.ttypes_to_use should be numeric'); end
  ttidx=P.ttypes_to_use;
  if any(ttidx<1 | ttidx>slength(M.ttype)), error('ttidx out of range'); end
  orig_npat = slength(M.pat); pidx = (1:orig_npat)'; M.pat.idx_orig = pidx;
  pidx = find(ismember(M.pat.ttype_idx,ttidx));
  M.pat = reorder_struct(M.pat,pidx); M.np = length(pidx);
  map = nan(orig_npat,1); map(M.pat.idx_orig) = (1:slength(M.pat));
  M.mut.pat_idx = map(M.mut.pat_idx);
  M.mut = reorder_struct_exclude(M.mut,isnan(M.mut.pat_idx));
  notes_fprintf('Using selected set of %d ttypes (%d patients).\n',length(ttidx),length(pidx));
end

% which patients to use
if ~isempty(P.patients_to_use)
  orig_npat = slength(M.pat); pidx = (1:orig_npat)'; M.pat.idx_orig = pidx;
  if isnumeric(P.patients_to_use)
    if ~isempty(P.ttype_grep)||~isempty(P.ttype_grepv) || ~isempty(P.ttypes_to_use)
      error('Please use string P.patients_to_use if using in conjunction with P.ttypes_to_use/ttype_grep/ttype_grepv');
    end
    pidx = P.patients_to_use;
  else
    pidx = find(ismember(M.pat.name,P.patients_to_use));
  end
  if any(pidx<1 | pidx>M.np), error('pidx out of range'); end
  M.pat = reorder_struct(M.pat,pidx); M.np = length(pidx);
  map = nan(orig_npat,1); map(M.pat.idx_orig) = (1:slength(M.pat));
  M.mut.pat_idx = map(M.mut.pat_idx);
  M.mut = reorder_struct_exclude(M.mut,isnan(M.mut.pat_idx));
  notes_fprintf('Using selected set of %d patients.\n',length(pidx));
end

if P.collapse_categories
  error('not yet implemented');
  % do nothing: will be handled later in code
end

% ignore indels?
if P.ignore_indels
  M.mut = reorder_struct_exclude(M.mut,M.mut.is_indel);
  notes_fprintf('Ignoring indels.\n');
end

% remove mutations that aren't assigned to a patient+gene
M.mut = reorder_struct_exclude(M.mut,isnan(M.mut.pat_idx)|isnan(M.mut.gene_idx));

if M.np==0 || slength(M.mut)==0
  fprintf('No data!\n');
  G=[];
  return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ng = M.ng;
np = M.np;
n_cov_models = length(M.cov.type);

if P.collapse_categories
  ncat = 1;
else
  ncat = M.ncat;
end

notes_fprintf('ng %d   np %d   ncat %d\n',ng,np,ncat);

% which genes to analyze
if ~exist('gta','var'), gta = 1:ng; end
if ischar(gta), gta={gta}; end
if iscellstr(gta), gta = listmap(gta,M.gene.name); gta(isnan(gta))=[]; end
if length(gta)<ng, notes_fprintf('  (analyzing %d genes)\n',length(gta)); end

% some speedups on M.mut
M.mut.is_nonsilent = (M.mut.is_coding & ~M.mut.is_silent);
M.mut = sort_struct(M.mut,'gene_idx');
[u ui uj] = unique(M.mut.gene_idx,'first');
h = histc(uj,1:length(u));
geneidx = cell(ng,1);
for i=1:length(u), geneidx{u(i)} = as_column(ui(i):ui(i)+h(i)-1);end

%%%%%%%%%%%%%%%%
% TERRITORY    %
%%%%%%%%%%%%%%%%

if P.increase_terr_to_max_cov_represented
  fprintf('Increasing territory to max coverage represented\n');
  if ~P.impute_full_cov_when_promotes_significance
    fprintf('\t(will have no effect on calculations, because P.impute_full_cov_when_promotes_significance = false)\n');
  end
  cov_idx = unique(M.pat.cov_idx);
  M.cov.gene_effect_terr2 = ceil(max([M.cov.gene_effect_terr M.cov.gene_effect_cov(:,cov_idx,:,:)],[],2));
else
  M.cov.gene_effect_terr2 = M.cov.gene_effect_terr;
end

%%%%%%%%%%%%%%%%
% TOTAL COUNTS %
%%%%%%%%%%%%%%%%

G=[];
G.gene = M.gene.name;
G.codelen = round(sum(sum(M.cov.gene_effect_terr2(:,1,:,2:5),4),3)/3);
flds = {'expr','exprmax','rt','RT','hiC','paz','nstrikes'};
for i=1:length(flds), if isfield(M.gene,flds{i}), G.(flds{i})=M.gene.(flds{i}); end, end

G.Nncd=0;
G.Nsyn=0;
G.Nmis=0;
G.Nnon=0;
G.Nspl=0;
for cov_idx=1:M.cov.ntype
  pidx = find(M.pat.cov_idx==cov_idx);
  G.Nncd = G.Nncd + sum(M.cov.gene_effect_cov(:,cov_idx,:,1),3)*length(pidx);
  G.Nsyn = G.Nsyn + sum(M.cov.gene_effect_cov(:,cov_idx,:,2),3)*length(pidx);
  G.Nmis = G.Nmis + sum(M.cov.gene_effect_cov(:,cov_idx,:,3),3)*length(pidx);
  G.Nnon = G.Nnon + sum(M.cov.gene_effect_cov(:,cov_idx,:,4),3)*length(pidx);
  G.Nspl = G.Nspl + sum(M.cov.gene_effect_cov(:,cov_idx,:,5),3)*length(pidx);
end
G.Nind = G.Nncd + G.Nsyn + G.Nmis + G.Nnon + G.Nspl;

G.Nncd = round(G.Nncd);
G.Nsyn = round(G.Nsyn);
G.Nmis = round(G.Nmis);
G.Nnon = round(G.Nnon);
G.Nspl = round(G.Nspl);
G.Nind = round(G.Nind);

G.nncd = as_column(histc(M.mut.gene_idx(strcmp('ncd',M.mut.effect)),1:M.ng));
G.nsyn = as_column(histc(M.mut.gene_idx(strcmp('syn',M.mut.effect)),1:M.ng));
G.nmis = as_column(histc(M.mut.gene_idx(strcmp('mis',M.mut.effect)),1:M.ng));
G.nnon = as_column(histc(M.mut.gene_idx(strcmp('non',M.mut.effect)),1:M.ng));
G.nspl = as_column(histc(M.mut.gene_idx(strcmp('spl',M.mut.effect)),1:M.ng));
G.nind = as_column(histc(M.mut.gene_idx(grepmi('indel_cod|indel_spl',M.mut.effect)),1:M.ng));

%%%%%%%%%%%%%%%
% TOTAL RATES %
%%%%%%%%%%%%%%%

globalrate_ncd = sum(G.nncd) / sum(G.Nncd);
globalrate_syn = sum(G.nsyn) / sum(G.Nsyn);
globalrate_mis = sum(G.nmis) / sum(G.Nmis);
globalrate_non = sum(G.nnon) / sum(G.Nnon);
globalrate_spl = sum(G.nspl) / sum(G.Nspl);
globalrate_ind = sum(G.nind) / sum(G.Nind);

notes_fprintf('Global rates (/Mb):  ncd %.2f   syn %.2f   mis %.2f  non %.2f  spl %.2f  ind %.2f\n',...
   globalrate_ncd*1e6,globalrate_syn*1e6,globalrate_mis*1e6,globalrate_non*1e6,globalrate_spl*1e6,globalrate_ind*1e6);

if globalrate_ncd==0 || globalrate_mis==0 || globalrate_syn==0 || globalrate_non==0 || globalrate_spl==0
  if P.scaling
    fprintf('WARNING: Disabling rate scaling because some globalrate(s)==0\n');
    P.scaling = false;
  end
end

%%%%%%%%%%%%%%%%%
% PATIENT RATES %
%%%%%%%%%%%%%%%%%

M.pat.Ntot = nan(np,1);
for c=1:M.cov.ntype
  pidx = find(M.pat.cov_idx==c);
  M.pat.Ntot(pidx) = fullsum(M.cov.gene_effect_cov(:,c,:,:));
end
M.pat.ntot = histc(M.mut.pat_idx,1:np);
M.pat.rtot = M.pat.ntot./M.pat.Ntot;
M.pat.log_rtot = max(-9,log10(M.pat.rtot));

if P.num_neighbor_patients>0
  % compute matrix of neighbor-patients
  if P.num_neighbor_patients>=20
    patneighb = zeros(np,np);
  else
    patneighb = sparse(np,np);
  end
  for p=1:np
    patdist = abs(M.pat.log_rtot-M.pat.log_rtot(p));
    [tmp ord] = sort(patdist);
    npidx = ord(1:min(length(ord),P.num_neighbor_patients+1));   % (include the patient itself)
    if isfield(M.pat,'callscheme')
      % restrict neighbor patients to within the same callscheme (avoids an N=0 bug)
      npidx = intersect(npidx,find(M.pat.callscheme==M.pat.callscheme(p)));
    end
    patneighb(p,npidx)=1;
  end
end

%%%%%%%%%%%%%%
% COVARIATES %
%%%%%%%%%%%%%%

fprintf('Processing covariates...\n');
nv = slength(M.V);
V = nan(ng,nv);
for vi=1:nv, V(:,vi) = M.V.val{vi}; end

% convert covariate raw values to Z-scores
Z = nan(ng,nv);
for vi=1:nv
  missing = isnan(V(:,vi)) | isinf(V(:,vi));
  mn = mean(V(~missing,vi));
  sd = std(V(~missing,vi),0);  % second parameter=0 means confirm default behavior of normalize by (N-1) not (N)
  Z(~missing,vi) = (V(~missing,vi)-mn)./sd;
end

if P.scaling
  notes_fprintf_once('Using scaling.\n');
else
  notes_fprintf_once('No scaling.\n');
end

notes_fprintf_once('still using original two-tailed hyge2cdf as stopping criterion\n');
notes_fprintf_once('for main bagels, using ncd+syn mutations\n');
notes_fprintf_once('for indel bagels, using all indels\n');

% will discover two sets of bagels:
%      main_bagel
%      indel_bagel

n_bageltypes = 2;
min_neighbors = [P.min_neighbors P.indel_min_neighbors];
max_neighbors = [P.max_neighbors P.indel_max_neighbors];
min_territory_ratio = [P.min_territory_ratio P.indel_min_territory_ratio];
qual_min = [P.qual_min P.indel_qual_min];

max_min_neighbors = max(min_neighbors);
max_max_neighbors = max(max_neighbors);
max_min_territory_ratio = max(min_territory_ratio);
min_qual_min = min(qual_min);

G.nnei = zeros(ng,1); G.bagel = nan(ng,max_max_neighbors);  % NOTE: not actuall used
G.nnei2 = G.nnei;
G = move_field_to_before(G,'nnei','nncd');
G = move_field_to_after(G,'nnei2','nnei');

notes_fprintf_once('Calculating p-value using 1D Projection method.\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROCESS EACH GENE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t1=0; t2=0; t3=0; t4=0;

G.npat = nan(ng,1);
G.nsite = nan(ng,1);
G.eff = nan(ng,1);
G.pmin = nan(ng,1);
G.pmid = nan(ng,1);
G.pmax = nan(ng,1);

if P.output_P0_table
  P0_table = nan(ng,np);
end


for g=as_row(gta), if ~mod(g,1000), fprintf('%d/%d ',g,ng); end

  tt = tic;

  gname = M.gene.name{g};

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % BAGEL FINDING
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
  % calculate distances from this gene
  df2 = bsxfun(@minus,Z,Z(g,:)).^2;
  dist2 = nansum(df2,2)./sum(~isnan(df2),2);
  [tmp,ord] = sort(dist2); ord = [g;ord(ord~=g)];

  per_type_nnei = zeros(n_bageltypes,1);
  per_type_bagel = repmat({nan(1,max_max_neighbors)},n_bageltypes,1);
  per_type_bagel_done = false(n_bageltypes,1);
  
  % expand bagel outward until quality falls below qual_min
  nfit=0; Nfit=0;
  for ni=0:max_max_neighbors, gidx = ord(ni+1);

    if P.scaling
      Ngene = G.Nncd(gidx)*(globalrate_ncd/globalrate_mis) + G.Nsyn(gidx)*(globalrate_syn/globalrate_mis);
      ngene = G.nncd(gidx) + G.nsyn(gidx);
    else
      Ngene = G.Nncd(gidx) + G.Nsyn(gidx);
      ngene = G.nncd(gidx) + G.nsyn(gidx);
    end

    if ni==0, ngene0=ngene; Ngene0=Ngene; end
    nfit=nfit+ngene; Nfit=Nfit+Ngene;
    
    % compare the gene being added to the central gene
    qual = 2*hyge2cdf(ngene,Ngene,ngene0,Ngene0);  % (two-sided)
    if qual>1, qual = 2-qual; end

    debug=DEBUG;
    if debug
      fprintf('%d %-10s  ngene %4d  Ngene %9d = %9.2d   ngene0 %4d  Ngene0 %9d = %9.2d    qual %9.2d\n',...
              ni,G.gene{gidx},ngene,round(Ngene),ngene/Ngene,ngene0,round(Ngene0),ngene0/Ngene0,qual);
    end

    % territory ratio
    territory_ratio = Nfit/Ngene0;
    
    % stopping criterion: stop if this gene would drop quality below qual_min
    if ((ni-1)>=max_min_neighbors && territory_ratio>=max_min_territory_ratio) && qual<min_qual_min, break; end
    
    % update gene's bagel
    if ni>0
      G.bagel(g,ni)=gidx;   % NOTE: not actually used anywhere
      % which categories include this bagel?
      for c=1:n_bageltypes
        if ~per_type_bagel_done(c) && (ni<=min_neighbors(c) || territory_ratio<min_territory_ratio(c) || (ni<=max_neighbors(c) && qual>=qual_min(c)))
          per_type_nnei(c) = per_type_nnei(c) + 1;
          per_type_bagel{c}(ni) = gidx;
        else
          per_type_bagel_done(c) = true;
        end
      end
    end      
  end % next neighborhood size

  % forced bagels?
  if ~isempty(P.forced_bagels)
    for i=1:length(P.forced_bagels)
      if strcmp(gname,P.forced_bagels{i}{1})
        fprintf('Using forced bagel for gene %s\n',gname);
        bagel = listmap(P.forced_bagels{i}{2},M.gene.name);
        bagel(isnan(bagel))=[];
        nb = length(bagel);
        per_type_bagel{1}(1:nb) = bagel; per_type_bagel{1}(nb+1:end)=nan;
        per_type_bagel{2} = per_type_bagel{1};
        per_type_nnei = [nb nb];
      end
    end      
  end

  G.nnei(g) = per_type_nnei(1);   % main bagel
  G.nnei2(g) = per_type_nnei(2);  % indel bagel

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % PROJECTION METHOD
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % STEP 0
  % retrieve counts for this gene + bagel

  % gene's mutation counts
  gene_midx = geneidx{g};
  midx1 = gene_midx(strcmp('mis',M.mut.effect(gene_midx)));
  gene_nmis = hist2d_fast(M.mut.pat_idx(midx1),M.mut.categ(midx1),1,np,1,ncat);
  midx2 = gene_midx(strcmp('non',M.mut.effect(gene_midx)));
  gene_nnon = hist2d_fast(M.mut.pat_idx(midx2),M.mut.categ(midx2),1,np,1,ncat);
  midx3 = gene_midx(strcmp('spl',M.mut.effect(gene_midx)));
  gene_nspl = hist2d_fast(M.mut.pat_idx(midx3),M.mut.categ(midx3),1,np,1,ncat);
  midx4 = gene_midx(grepm('indel_(cod|spl)',M.mut.effect(gene_midx)));
  if ~isempty(midx4), gene_nind = as_column(histc(M.mut.pat_idx(midx4),1:np)); else gene_nind = zeros(np,1); end

  % find # unique patients and sites (only for display in table; not used directly in calculation)
  midx = [midx1;midx2;midx3;midx4];
  if isempty(midx)
    G.nsite(g) = 0;
    G.npat(g) = 0;
  else
    G.nsite(g) = length(unique_combos(M.mut.chr(midx),M.mut.start(midx)));
    G.npat(g) = length(unique(M.mut.patient(midx)));
  end

  % gene's coverage (or territory)
  if P.impute_full_cov_when_promotes_significance
    % use territory
    tmp = squeeze(M.cov.gene_effect_terr2(g,1,:,2:5));
    gene_Nsyn = repmat(tmp(:,1)',np,1);
    gene_Nmis = repmat(tmp(:,2)',np,1);
    gene_Nnon = repmat(tmp(:,3)',np,1);
    gene_Nspl = repmat(tmp(:,4)',np,1);
  else  % use coverage
    gene_Nsyn = zeros(np,ncat);
    gene_Nmis = zeros(np,ncat);
    gene_Nnon = zeros(np,ncat);
    gene_Nspl = zeros(np,ncat);
    for c=1:ncat
      gene_Nsyn(:,c) = squeeze(M.cov.gene_effect_cov(g,M.pat.cov_idx,c,2));
      gene_Nmis(:,c) = squeeze(M.cov.gene_effect_cov(g,M.pat.cov_idx,c,3));
      gene_Nnon(:,c) = squeeze(M.cov.gene_effect_cov(g,M.pat.cov_idx,c,4));
      gene_Nspl(:,c) = squeeze(M.cov.gene_effect_cov(g,M.pat.cov_idx,c,5));
    end
  end
  gene_Nind = sum(gene_Nsyn + gene_Nmis + gene_Nnon + gene_Nspl,2);

  % main sphere (=gene+bagel): noncoding and synonymous
  main_sphere_nncd = zeros(np,ncat);
  main_sphere_Nncd = zeros(np,ncat);
  main_sphere_nsyn = zeros(np,ncat);
  main_sphere_Nsyn = zeros(np,ncat);
  bageltype = 1;
  main_sphere = [g per_type_bagel{bageltype}(1:per_type_nnei(bageltype))];
  main_sphere_midx = cat(1,geneidx{main_sphere});
  for c=1:ncat
    midx = main_sphere_midx(strcmp('ncd',M.mut.effect(main_sphere_midx)) & M.mut.categ(main_sphere_midx)==c);
    if ~isempty(midx), main_sphere_nncd(:,c) = as_column(histc(M.mut.pat_idx(midx),1:np)); end
    midx = main_sphere_midx(strcmp('syn',M.mut.effect(main_sphere_midx)) & M.mut.categ(main_sphere_midx)==c);
    if ~isempty(midx), main_sphere_nsyn(:,c) = as_column(histc(M.mut.pat_idx(midx),1:np)); end
    for cov_idx=1:M.cov.ntype
      pidx = (M.pat.cov_idx==cov_idx);
      main_sphere_Nncd(pidx,c) = sum(M.cov.gene_effect_cov(main_sphere,cov_idx,c,1),1);
      main_sphere_Nsyn(pidx,c) = sum(M.cov.gene_effect_cov(main_sphere,cov_idx,c,2),1);
    end
  end
  % expand to neighboring patients (if specified)
  if P.num_neighbor_patients>0
    main_sphere_nncd = patneighb*main_sphere_nncd;
    main_sphere_Nncd = patneighb*main_sphere_Nncd;
    main_sphere_nsyn = patneighb*main_sphere_nsyn;
    main_sphere_Nsyn = patneighb*main_sphere_Nsyn;
  end

  % indel bagel
  indel_bagel_nind = zeros(np,1);
  indel_bagel_Nind = zeros(np,1);
  bageltype = 2;
  indel_bagel = per_type_bagel{bageltype}(1:per_type_nnei(bageltype));
  if isempty(indel_bagel), error('what?'); end
  indel_bagel_midx = cat(1,geneidx{indel_bagel});
  midx = indel_bagel_midx(strncmp('indel',M.mut.effect(indel_bagel_midx),5));
  if ~isempty(midx), indel_bagel_nind = as_column(histc(M.mut.pat_idx(midx),1:np)); end
  for cov_idx=1:M.cov.ntype
    pidx = (M.pat.cov_idx==cov_idx);
    indel_bagel_Nind(pidx) = sum(sum(sum(M.cov.gene_effect_cov(indel_bagel,cov_idx,:,:),3),4),1);
  end

  % STEP 1
  % compute probability of seeing zero mutations in each degree

  if P.scaling
    main_sphere_Nncd = main_sphere_Nncd * (globalrate_ncd/globalrate_mis);
    main_sphere_Nsyn = main_sphere_Nsyn * (globalrate_syn/globalrate_mis);
    gene_Nnon = gene_Nnon * (globalrate_non/globalrate_mis);
    gene_Nspl = gene_Nspl * (globalrate_spl/globalrate_mis);
  end

  main_ntot = gene_nmis + gene_nnon + gene_nspl + main_sphere_nncd + main_sphere_nsyn;
  main_Ntot = gene_Nmis + gene_Nnon + gene_Nspl + main_sphere_Nncd + main_sphere_Nsyn;
  indel_ntot = gene_nind + indel_bagel_nind;
  indel_Ntot = gene_Nind + indel_bagel_Nind;

  t1=t1+toc(tt); tt = tic;

              %  1   2   3   4
  ndeg = 4;   % mis non spl ind
  P0 = nan(np,ndeg);
  has_mutation = false(np,ndeg);
  for d=1:ndeg
    if d<4
      n_total = main_ntot; N_total = main_Ntot;
      if d==1
        n_signal = gene_nmis; N_signal = gene_Nmis;
      elseif d==2
        n_signal = gene_nnon; N_signal = gene_Nnon;
      else % d==3
        n_signal = gene_nspl; N_signal = gene_Nspl;
      end
    else % d==4
      n_total = indel_ntot; N_total = indel_Ntot;
      n_signal = gene_nind; N_signal = gene_Nind;
    end
    % make sure we never have a more than 1000x difference between Nsignal and Ntotal: is not realistic
    N_total(N_total>1000*N_signal) = 1000*N_signal(N_total>1000*N_signal);     % to catch an N=0 bug
    % make sure values are rounded
    N_total = round(N_total); n_total = round(n_total); N_signal = round(N_signal);
    % make sure we never have N>2*n
    N_signal(2*n_signal>N_signal)=2*n_signal(2*n_signal>N_signal);
    N_total(2*n_total>N_total)=2*n_total(2*n_total>N_total);
    % find probabilities
    t2=t2+toc(tt); tt = tic;
    p = my_hygepdf(0,N_total,n_total,N_signal);  % (2D vectorization was being screwed up by all the error checking in hygepdf.m)
    t3=t3+toc(tt); tt = tic;
    P0(:,d) = prod(p,2);
    has_mutation(:,d) = sum(n_signal,2)>0;
    t2=t2+toc(tt); tt = tic;
  end
  P0(isnan(P0))=1;
  P0(P0>1)=1;
  P0(P0<0)=0;
  P1 = 1-P0;

  t2=t2+toc(tt); tt = tic;

  % for each sample, prioritize mutation categories according to how likely
  % it would be for this gene x sample to have a mutation there by chance.
  %
  % determine each patient's priority order of categories (according to P1)
  %  left column of "priority" =  lowest priority =  most likely category of mutation
  % right column of "priority" = highest priority = least likely category of mutation
  [tmp priority] = sort(P1,2,'descend');
  % sort the P arrays to put the columns in least->most priority order
  shft = (priority - repmat(1:ndeg,np,1));
  map = reshape(1:(np*ndeg),np,ndeg);
  newmap = map + shft*np;
  P0 = P0(newmap);
  P1 = P1(newmap);
  has_mutation = has_mutation(newmap);

  % STEP 2
  % for each sample, compute probability that it would have been of each degree.
  % where d=0 (no mut) ..... ndeg (most extreme mut)
  Pdeg = nan(np,ndeg+1);
  for d=0:ndeg
    % has to be clear in all categories > d
    Pdeg(:,d+1) = prod(P0(:,d+1:end),2);
    % and nonclear in category d (if d>0)
    if d>0, Pdeg(:,d+1) = Pdeg(:,d+1) .* P1(:,d); end
  end

  if P.output_P0_table
    P0_table(g,:) = Pdeg(:,1);
  end

  %% STEP 2a: calculate score for a sample being of each possible degree
  %% (uses new style, where score = -log10 probability of having a mutation in that category
  %% (zero score for no mutations)
  Sdeg = [zeros(np,1) -log10(P1)];

  % null score boost
  priority2 = [zeros(np,1) priority];
  Sdeg(priority2>=2) = Sdeg(priority2>=2) + P.null_score_boost;

  % score multiplier (to increase resolution of scores)
  Sdeg = Sdeg * P.score_multiplier;
  notes_fprintf_once('using P.score_multiplier = %d\n',P.score_multiplier);

  % round all scores UP to nearest integer
  Sdeg = ceil(Sdeg);
  notes_fprintf_once('rounding all scores UP to nearest integer\n');

  % score cap
  Sdeg = min(Sdeg,1000);   % because inf will crash the C version!
  
  % STEP 3
  % determine actual degree of each sample, and score_obs for the gene

  degree = zeros(np,1);
  score_obs = 0;
  for p = 1:np
    for d = ndeg:-1:1
      if has_mutation(p,d), degree(p) = d; break; end
    end
    score_obs = score_obs + Sdeg(p,degree(p)+1);
  end

  % impose minimum effect size by decreasing score_obs
  score_obs_unreduced = score_obs;
  if length(P.min_effect_size)==ng
    score_obs = score_obs / P.min_effect_size(g);
  elseif length(P.min_effect_size)==1
    score_obs = score_obs / P.min_effect_size;
  else
    error('invalid P.min_effect_size');
  end

  % for zero score, don't bother doing convolutions
%  if score_obs<=0, G.p(g)=1; continue; end
%  notes_fprintf_once('skipping genes with zero mutations!\n');
  notes_fprintf_once('no longer skipping genes with zero mutations!\n');

  % STEP 4
  % compute P value for gene by convolutions
  % --> 95% of the total runtime is in this step (t4/ttot)

  if P.numbins_scheme==1
    numbins = max(10,score_obs*2);
  elseif P.numbins_scheme==2
    numbins = ceil(score_obs+max(5,0.2*score_obs));
  else
    error('unknown P.numbins_scheme');
  end

  % allocate space for convolutions
  try
    H = zeros(numbins,1);
    newH = zeros(numbins,ndeg+1);
  catch me
    fprintf('WEIRD ERROR!>>>>>'); keyboard;
  end

  % compute score_exp and effect_size
  score_exp = sum(Sdeg(:).*Pdeg(:));
  if score_obs_unreduced==0 && score_exp==0
    effect_size = 0;
  else
    effect_size = score_obs_unreduced/score_exp;
  end
  G.eff(g) = effect_size;

  [pmax pmin] = projection_1d_convolutions_fast(Sdeg,Pdeg,score_obs,numbins,H,newH);
  % pmax is the p-value generated by the standard procedure we've always used.
  % pmin is the (better) p-value that goes one discrete step further
  % to solve the problem of discrete statistics, we want to take a randomly chosen position between pmin and pmax.
  G.pmax(g) = pmax;        % for sorting genelist
  G.pmid(g) = pmin+rand*(pmax-pmin);  % for q-q plot
  G.pmin(g) = pmin;        % for future reference

  t4=t4+toc(tt); tt = tic;
    
  ttot = t1+t2+t3+t4;
    
  progress_fprintf('%5d %12s %4.1f  %12.10f %12.10f %12.10f   %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n',...
          g,M.gene.name{g},G.eff(g),G.pmin(g),G.pmid(g),G.pmax(g),t1,t2,t3,t4,ttot,ttot/g);

  debug=DEBUG;
  if debug
    has_idx=find(any(has_mutation,2));
    hasnt_idx = find(~any(has_mutation,2));
    pidx = [has_idx(1:min(length(has_idx),13));hasnt_idx(1:min(length(hasnt_idx),2))];
    x=[];x.patient=M.pat.name;x.callscheme=M.pat.callscheme_name;
    x.main_ntot=main_ntot;x.indel_ntot=indel_ntot;x.main_Ntot=round(main_Ntot);
    x.priority=priority;x.has_mutation=has_mutation;x.degree=degree;x.P0=P0;x.Pdeg=Pdeg;x.Sdeg=Sdeg;
    pr(x,pidx)
    keyboard
  end
  
end, fprintf('\n');   % next gene

% FDR
fprintf('%s [sort]\n',datestr(now));
G = sort_struct(G,'pmax');
fprintf('%s [FDR]\n',datestr(now));
G.q = calc_fdr_value(G.pmax);

% remove some columns
G = rmfield_if_exist(G,{'bagel'});  % too inconvenient for display
G = rmfield(G,grep('^N',fieldnames(G)));

% SAVE+RETURN RESULTS
if exist('outdir','var'), save([outdir '/results.mat'],'G'); end
if exist('outdir','var'), save_struct(G,[outdir '/sig_genes.txt']); end
pr(G,1:min(30,length(gta)))

if P.output_P0_table
  filename = [outdir '/P0_table.mat'];
  fprintf('Writing P0 table to %s\n',filename);
  P0_table_genes = M.gene.name;
  P0_table_patients = M.pat.name;
  save(filename,'P0_table','P0_table_genes','P0_table_patients');
end

nsg = sum(G.q<=0.1);
notes_fprintf('%d genes with q<=0.1',nsg);
idx = grep('^OR\d',G.gene(G.q<0.3),1);
if ~isempty(idx)
  notes_fprintf(', with OR');
  if length(idx)>1, notes_fprintf('s'); end
  notes_fprintf(' at ');
  for i=1:length(idx), notes_fprintf('%d',idx(i)); if i<length(idx), notes_fprintf(', '); end, end
end
notes_fprintf('\n');

if exist('progress_f','var'), fclose(progress_f); end
if exist('notes_f','var'), fclose(notes_f); end

fprintf('%s [bottom]\n',datestr(now));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  function progress_fprintf(varargin)
    fprintf(varargin{:});
    if exist('progress_f','var')
      fmt = regexprep(varargin{1},'(\S)\s+(\S)','$1\t$2');  % replace spaces with tabs
      fprintf(progress_f,fmt,varargin{2:end});
    end
  end
 
  function notes_fprintf(varargin)
    fprintf(varargin{:});
    if exist('notes_f','var'), fprintf(notes_f,varargin{:}); end
  end

  function notes_fprintf_once(varargin)
    hashval = sum(double(varargin{1}));
    if length(flaghash)>=hashval && flaghash(hashval)==1, return; end
    notes_fprintf(varargin{:});
    flaghash(hashval)=1;
  end

end
