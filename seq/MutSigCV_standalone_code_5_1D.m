function G = MutSigCV_standalone_code_5_1D(M,gta,P,outdir)

% using COMPACT representations

% with 1D projection

% RESTORE THE USE OF BAGELS

flaghash = sparse([]);

% omitted gta?
if nargin==2 && (isempty(gta)||isstruct(gta)), P=gta; clear gta; end
if nargin==3 && (isempty(gta)||isstruct(gta)), outdir=P; P=gta; clear gta; end

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'collapse_categories',false);
P = impose_default_value(P,'ttypes_to_use',[]);
P = impose_default_value(P,'patients_to_use',[]);
P = impose_default_value(P,'categories_to_use',[]);
P = impose_default_value(P,'ignore_indels',false);
P = impose_default_value(P,'demote_nulls',false);
P = impose_default_value(P,'special_treatment_for_category',[]);
P = impose_default_value(P,'use_bagels',true);
P = impose_default_value(P,'use_bagels_from',[]);
P = impose_default_value(P,'per_category_bagels',false);
P = impose_default_value(P,'min_neighbors',0);
P = impose_default_value(P,'max_neighbors',50);
P = impose_default_value(P,'qual_min',0.05);
P = impose_default_value(P,'numbins_scheme',2);
P = impose_default_value(P,'use_flanking_mutations',false);
P = impose_default_value(P,'use_silent_mutations',true);
if (~P.use_silent_mutations), error('not implemented'); end
P = impose_default_value(P,'scale_flanking_mutations',true);
P = impose_default_value(P,'scale_silent_mutations',false);
P = impose_default_value(P,'signal','nonsilent');
P = impose_default_value(P,'score_multiplier',10);

if P.use_flanking_mutations
  P = impose_default_value(P,'total','silent+nonsilent+flanking+bagel_silent+bagel_flanking');
else
  P = impose_default_value(P,'total','silent+nonsilent+bagel_silent');
end

if exist('outdir','var')
  ede(outdir);
  progress_f = fopen([outdir '/progress.txt'],'wt');
  notes_f = fopen([outdir '/notes.txt'],'wt');
end

if isfield(M,'pat')&&~isfield(M,'patient'),M=rename_field(M,'pat','patient');end

demand_fields(M,{'mut','patient','gene','cov','ng','np'});
demand_fields(M.patient,{'cov_idx'});
demand_fields(M.mut,{'gene_idx','pat_idx','is_coding','is_flank','is_silent'});
demand_fields(M.cov,{'gene_flk_cov','gene_sil_cov','gene_non_cov'});
demand_fields(M.cov,{'gene_flk_cov_tot','gene_sil_cov_tot','gene_non_cov_tot'});

if ~isempty(P.ttypes_to_use) && ~isempty(P.patients_to_use)
  error('please use only one');
end

% which ttypes to use
if ~isempty(P.ttypes_to_use)
  if ~isnumeric(P.ttypes_to_use), error('P.ttypes_to_use should be numeric'); end
  ttidx=P.ttypes_to_use;
  if any(ttidx<1 | ttidx>slength(M.ttype)), error('ttidx out of range'); end
  pidx = find(ismember(M.patient.ttype_idx,ttidx));
  M.patient = reorder_struct(M.patient,pidx); M.np = length(pidx);
  M.mut = reorder_struct(M.mut,ismember(M.mut.pat_idx,pidx));
  M.mut.pat_idx = listmap(M.mut.patient,M.patient.name);
  notes_fprintf('Using selected set of %d ttypes (%d patients).\n',length(ttidx),length(pidx));
end

% which patients to use
if ~isempty(P.patients_to_use)
  if ~isnumeric(P.patients_to_use), error('P.patients_to_use should be numeric'); end
  pidx=P.patients_to_use;
  if any(pidx<1 | pidx>M.np), error('pidx out of range'); end
  M.patient = reorder_struct(M.patient,pidx); M.np = length(pidx);
  M.mut = reorder_struct(M.mut,ismember(M.mut.pat_idx,pidx));
  M.mut.pat_idx = listmap(M.mut.patient,M.patient.name);
  notes_fprintf('Using selected set of %d patients.\n',length(pidx));
end

% which categories to analyze
if ~isempty(P.categories_to_use) && P.collapse_categories
  error('inconsistent P settings');
end

if P.collapse_categories
  % do nothing: will be handled later in code
end

if ~isempty(P.categories_to_use)
  if ~isnumeric(P.categories_to_use), error('P.categories_to_use should be numeric'); end
  cidx=P.categories_to_use;
  if any(cidx<1 | cidx>M.cov.ncat), error('cidx out of range'); end
  f=grep('^gene_.*_terr$',fieldnames(M.cov));for i=1:length(f),M.cov.(f{i})=M.cov.(f{i})(:,cidx);end
  f=grep('^gene_.*_cov$',fieldnames(M.cov));for i=1:length(f),M.cov.(f{i})=M.cov.(f{i})(:,:,cidx);end
  M.cov.ncat = length(cidx); 
  M.cov.categ=reorder_struct(M.cov.categ,cidx);
  M.ncat=M.cov.ncat; M.categ = M.cov.categ;
  M.mut.categ = mapacross(M.mut.categ,cidx,1:length(cidx));
  M.mut.categ_ignoring_null_categ = mapacross(M.mut.categ_ignoring_null_categ,cidx,1:length(cidx));
  notes_fprintf('Using selected set of %d categories.\n',length(cidx));
end

% is there a null category?
if ~isempty(grepi('null',M.categ.name))
  notes_fprintf('Using a null category.\n');
else
  notes_fprintf('Not using a null category.\n');
  P.demote_nulls = true;
end

if P.demote_nulls
  notes_fprintf('Demoting nulls.\n');
  M.mut.categ = M.mut.categ_ignoring_null_categ;
end

% ignore indels?
indel_idx = grep('INS|DEL|Indel',M.mut.classification,1);
if P.ignore_indels
  M.mut = reorder_struct_exclude(M.mut,indel_idx);
  notes_fprintf('Ignoring indels.\n');
else
  if ~isempty(indel_idx)
    if ~isempty(grepi('indel',M.categ.name)) || M.ncat==1
      notes_fprintf('Keeping indels.\n');
    else
      error('Don''t know what category to put indels in');
    end
  end
end

% remove mutations that aren't assigned to a patient+gene+category
M.mut = reorder_struct_exclude(M.mut,isnan(M.mut.pat_idx)|isnan(M.mut.gene_idx)|isnan(M.mut.categ));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ng = M.ng;
np = M.np;
n_cov_models = size(M.cov.gene_non_cov,2);

if P.collapse_categories
  ncat = 1;
else
  ncat = M.ncat;
end
ncols = ncat+1;

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
% TOTAL COUNTS %
%%%%%%%%%%%%%%%%

G=[];
G.gene = M.gene.name;

G.Nnon=0;
G.Nsil=0;
G.Nflk=0;
for cov_idx=1:n_cov_models
  pidx = find(M.patient.cov_idx==cov_idx);
  G.Nnon = G.Nnon + M.cov.gene_non_cov_tot(:,cov_idx)*length(pidx);
  G.Nsil = G.Nsil + M.cov.gene_sil_cov_tot(:,cov_idx)*length(pidx);
  G.Nflk = G.Nflk + M.cov.gene_flk_cov_tot(:,cov_idx)*length(pidx);
end
G.Nnon = round(G.Nnon);
G.Nsil = round(G.Nsil);
G.Nflk = round(G.Nflk);

G.nnon = histc(M.mut.gene_idx(M.mut.is_nonsilent),1:M.ng);
G.nsil = histc(M.mut.gene_idx(M.mut.is_silent),1:M.ng);
G.nflk = histc(M.mut.gene_idx(M.mut.is_flank),1:M.ng);

%%%%%%%%%%%%%%%
% TOTAL RATES %
%%%%%%%%%%%%%%%

globalrate_non = sum(G.nnon) / sum(G.Nnon);
globalrate_sil = sum(G.nsil) / sum(G.Nsil);
globalrate_flk = sum(G.nflk) / sum(G.Nflk);

notes_fprintf('Global rates (/Mb):  non %.2f   sil %.2f   flank %.2f\n',...
   globalrate_non*1e6,globalrate_sil*1e6,globalrate_flk*1e6);

%%%%%%%%%%%%%%
% BAGELS     %
%%%%%%%%%%%%%%

if ~P.use_bagels && ~isempty(P.use_bagels_from)
  error('incompatible parameters');

elseif ~P.use_bagels
  notes_fprintf('Not using bagels.\n');

elseif P.use_bagels && ~isempty(P.use_bagels_from)
  fprintf('Loading existing bagels from %s\n',P.use_bagels_from);
  tmp = load(P.use_bagels_from);
  G.bagel = tmp.G.bagel;
  G.nnei = tmp.G.nnei;

elseif P.use_bagels
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
  
  min_neighbors = P.min_neighbors;
  max_neighbors = P.max_neighbors;
  qual_min = P.qual_min;

  if P.per_category_bagels
    notes_fprintf('   using per-category bagels\n');
    if length(min_neighbors)==1, min_neighbors = repmat(min_neighbors,ncat,1); end
    if length(max_neighbors)==1, max_neighbors = repmat(max_neighbors,ncat,1); end
    if length(qual_min)==1, qual_min = repmat(qual_min,ncat,1); end
    if length(min_neighbors)~=ncat || length(max_neighbors)~=ncat || length(qual_min)~=ncat
      error('mismatched parameter lengths');
    end
    M.categ.min_neighbors = as_column(min_neighbors);
    M.categ.max_neighbors = as_column(max_neighbors);
    M.categ.qual_min = as_column(qual_min);
    notes_fprintf('%s',pr(M.categ));
  else
    if length(min_neighbors)>1 || length(max_neighbors)>1 || length(qual_min)>1
      error('please set P.per_category_bagels=true');
    end
    notes_fprintf('   min_neighbors %d   max_neighbors %d   qual_min %d\n',min_neighbors,max_neighbors,qual_min);
  end

  notes_fprintf_once('still using original two-tailed hyge2cdf as stopping criterion\n');
  notes_fprintf_once('for bagels, using:  ');
  if P.use_flanking_mutations
    notes_fprintf_once('flanking');
    if P.scale_flanking_mutations
      notes_fprintf_once(' (scaled)');
    else
      notes_fprintf_once(' (unscaled)');
    end
  end
  if P.use_silent_mutations
    if P.use_flanking_mutations
      notes_fprintf_once(' + ');
    end
    notes_fprintf_once('silent');
    if P.scale_silent_mutations
      notes_fprintf_once(' (scaled)');
    else
      notes_fprintf_once(' (unscaled)');
    end
  end

  max_min_neighbors = max(min_neighbors);
  max_max_neighbors = max(max_neighbors);
  min_qual_min = min(qual_min);

  G.nnei = zeros(ng,1); G.bagel = nan(ng,max_max_neighbors);
  G.Nsil_bagel = zeros(ng,1); G.Nflk_bagel = zeros(ng,1);
  G.nsil_bagel = zeros(ng,1); G.nflk_bagel = zeros(ng,1);
  if P.per_category_bagels
    per_cat_nnei = zeros(ng,ncat);
    per_cat_bagel = repmat({nan(ng,max_max_neighbors)},ncat,1);
    per_cat_bagel_done = false(ng,ncat);
  end

  fprintf('Finding bagels...\n');
  for g=as_row(gta), if ~mod(g,1000), fprintf('%d/%d ',g,ng); end
    
    % calculate distances from this gene
    df2 = bsxfun(@minus,Z,Z(g,:)).^2;
    dist2 = nansum(df2,2)./sum(~isnan(df2),2);
    [tmp,ord] = sort(dist2); ord = [g;ord(ord~=g)];
    
    % expand bagel outward until quality falls below qual_min
    nfit=0; Nfit=0;
    for ni=0:max_max_neighbors, gidx = ord(ni+1);

      Ngene = 0;
      ngene =0;
      if P.use_flanking_mutations
        ngene = ngene + G.nflk(gidx);
        if P.scale_flanking_mutations
          Ngene = Ngene + round(G.Nflk(gidx)/(globalrate_non/globalrate_flk));
        else
          Ngene = Ngene + G.Nflk(gidx);
        end
      end
      if P.use_silent_mutations
        ngene = ngene + G.nsil(gidx);
        if P.scale_silent_mutations
          Ngene = Ngene + round(G.Nsil(gidx)/(globalrate_non/globalrate_sil));
        else
          Ngene = Ngene + G.Nsil(gidx);
        end
      end
      
      if ni==0, ngene0=ngene; Ngene0=Ngene; end
      nfit=nfit+ngene; Nfit=Nfit+Ngene;
      
      % compare the gene being added to the central gene
      qual = 2*hyge2cdf(ngene,Ngene,ngene0,Ngene0);  % (two-sided)
      if qual>1, qual = 2-qual; end

      debug=true;
      if debug
        fprintf('%d %-10s  ngene %4d  Ngene %9d = %9.2d   ngene0 %4d  Ngene0 %9d = %9.2d    qual %9.2d\n',...
                ni,G.gene{gidx},ngene,round(Ngene),ngene/Ngene,ngene0,round(Ngene0),ngene0/Ngene0,qual);
      end
        
      % stopping criterion: stop if this gene would drop quality below qual_min
      if G.nnei(g)>=max_min_neighbors && qual<min_qual_min, break; end
      
      % update gene's bagel
      G.nnei(g) = ni;
      if ni>0
        G.bagel(g,ni)=gidx;
        % compute totals to display in G table (for display only)
        G.Nsil_bagel(g) = G.Nsil_bagel(g) + G.Nsil(gidx);
        G.Nflk_bagel(g) = G.Nflk_bagel(g) + G.Nflk(gidx);
        G.nsil_bagel(g) = G.nsil_bagel(g) + G.nsil(gidx);
        G.nflk_bagel(g) = G.nflk_bagel(g) + G.nflk(gidx);
        if P.per_category_bagels
          % which categories include this bagel?
          for c=1:ncat
            if ~per_cat_bagel_done(g,c) && (ni<=min_neighbors(c) || (ni<=max_neighbors(c) && qual>=qual_min(c)))
              per_cat_nnei(g,c) = per_cat_nnei(g,c) + 1;
              per_cat_bagel{c}(g,ni) = gidx;
            else
              per_cat_bagel_done(g,c) = true;
            end  
          end
        end
      end      
    end % next neighborhood size
  end, fprintf('\n'); % next gene
  
end

%%%%%%%%%%%%%%
% PROJECTION %
%%%%%%%%%%%%%%

signal = P.signal;
total = P.total;

notes_fprintf('Calculating p-value using 1D Projection method...\n');

t1=0; t2=0; t3=0; t4=0;

G.pmin = nan(ng,1);
G.p = nan(ng,1);
G.pmax = nan(ng,1);

for g=as_row(gta)

  tt = tic;

  % STEP 0
  % retrieve counts for this gene
  % nsignal, Nsignal   (np,ncat)   typically this is nonsilent of the gene itself
  % ntotal,  Ntotal    (np,ncat)   e.g. nonsilent+silent+flanking of the gene itself, plus silent+flanking of the bagel

  gene_midx = geneidx{g};

  % gene nonsilent
  midx = gene_midx(M.mut.is_nonsilent(gene_midx));
  if P.collapse_categories
    if isempty(midx)
      gene_nnon = zeros(np,1);
    else
      gene_nnon = as_column(histc(M.mut.pat_idx(midx),1:np));
    end
    gene_Nnon = as_column(M.cov.gene_non_cov_tot(g,M.patient.cov_idx));
  else
    gene_nnon = hist2d_fast(M.mut.pat_idx(midx),M.mut.categ(midx),1,np,1,ncat);
    gene_Nnon = zeros(np,ncat); for c=1:ncat, gene_Nnon(:,c) = M.cov.gene_non_cov(g,M.patient.cov_idx,c)'; end
  end
  % gene silent
  midx = gene_midx(M.mut.is_silent(gene_midx));
  if P.collapse_categories
    if isempty(midx)
      gene_nsil = zeros(np,1);
    else
      gene_nsil = as_column(histc(M.mut.pat_idx(midx),1:np));
    end
    gene_Nsil = as_column(M.cov.gene_sil_cov_tot(g,M.patient.cov_idx));
  else
    gene_nsil = hist2d_fast(M.mut.pat_idx(midx),M.mut.categ(midx),1,np,1,ncat);
    gene_Nsil = zeros(np,ncat); for c=1:ncat, gene_Nsil(:,c) = M.cov.gene_sil_cov(g,M.patient.cov_idx,c)'; end
  end
  % gene flanking
  midx = gene_midx(M.mut.is_flank(gene_midx));
  if P.collapse_categories
    if isempty(midx)
      gene_nflk = zeros(np,1);
    else
      gene_nflk = as_column(histc(M.mut.pat_idx(midx),1:np));
    end
    gene_Nflk = as_column(M.cov.gene_flk_cov_tot(g,M.patient.cov_idx));
  else
    gene_nflk = hist2d_fast(M.mut.pat_idx(midx),M.mut.categ(midx),1,np,1,ncat);
    gene_Nflk = zeros(np,ncat); for c=1:ncat, gene_Nflk(:,c) = M.cov.gene_flk_cov(g,M.patient.cov_idx,c)'; end
  end

  if P.use_bagels
    if ~P.per_category_bagels
      bagel = G.bagel(g,1:G.nnei(g));
      bagel_midx = cat(1,geneidx{bagel});
    else % P.per_category_bagels
      for c=1:ncat
        bagel{c} = per_cat_bagel{c}(g,1:per_cat_nnei(g,c));
        bagel_midx{c} = cat(1,geneidx{bagel{c}});
      end
    end
  else
    bagel = [];
  end

  if ~isempty(bagel) 
    % bagel nonsilent
    if P.collapse_categories
      midx = bagel_midx(M.mut.is_nonsilent(bagel_midx));
      if isempty(midx)
        bagel_nnon = zeros(np,1);
      else
        bagel_nnon = as_column(histc(M.mut.pat_idx(midx),1:np));
      end
      bagel_Nnon = zeros(np,1);
      for cov_idx=1:n_cov_models
        pidx = (M.patient.cov_idx==cov_idx);
        bagel_Nnon(pidx) = sum(M.cov.gene_non_cov_tot(bagel,cov_idx),1);
      end
    elseif ~P.per_category_bagels
      midx = bagel_midx(M.mut.is_nonsilent(bagel_midx));
      bagel_nnon = hist2d_fast(M.mut.pat_idx(midx),M.mut.categ(midx),1,np,1,ncat);
      bagel_Nnon = zeros(np,ncat);
      for cov_idx=1:n_cov_models
        pidx = (M.patient.cov_idx==cov_idx);
        bagel_Nnon(pidx,:) = repmat(shiftdim(sum(M.cov.gene_non_cov(bagel,cov_idx,:),1))',sum(pidx),1);
      end
    else % P.per_category_bagels
      bagel_nnon = zeros(np,ncat);
      bagel_Nnon = zeros(np,ncat);
      for c=1:ncat
        midx = bagel_midx{c}(M.mut.is_nonsilent(bagel_midx{c}) & M.mut.categ(bagel_midx{c})==c);
        if ~isempty(midx)
          bagel_nnon(:,c) = as_column(histc(M.mut.pat_idx(midx),1:np));
        end
        for cov_idx=1:n_cov_models
          pidx = (M.patient.cov_idx==cov_idx);
          bagel_Nnon(pidx,c) = sum(M.cov.gene_non_cov(bagel{c},cov_idx,c),1);
        end
      end
    end

    % bagel silent
    if P.collapse_categories
      midx = bagel_midx(M.mut.is_silent(bagel_midx));
      if isempty(midx)
        bagel_nsil = zeros(np,1);
      else
        bagel_nsil = as_column(histc(M.mut.pat_idx(midx),1:np));
      end
      bagel_Nsil = zeros(np,1);
      for cov_idx=1:n_cov_models
        pidx = (M.patient.cov_idx==cov_idx);
        bagel_Nsil(pidx) = sum(M.cov.gene_sil_cov_tot(bagel,cov_idx),1);
      end
    elseif ~P.per_category_bagels
      midx = bagel_midx(M.mut.is_silent(bagel_midx));
      bagel_nsil = hist2d_fast(M.mut.pat_idx(midx),M.mut.categ(midx),1,np,1,ncat);
      bagel_Nsil = zeros(np,ncat);
      for cov_idx=1:n_cov_models
        pidx = (M.patient.cov_idx==cov_idx);
        bagel_Nsil(pidx,:) = repmat(shiftdim(sum(M.cov.gene_sil_cov(bagel,cov_idx,:),1))',sum(pidx),1);
      end
    else % P.per_category_bagels
      bagel_nsil = zeros(np,ncat);
      bagel_Nsil = zeros(np,ncat);
      for c=1:ncat
        midx = bagel_midx{c}(M.mut.is_silent(bagel_midx{c}) & M.mut.categ(bagel_midx{c})==c);
        if ~isempty(midx)
          bagel_nsil(:,c) = as_column(histc(M.mut.pat_idx(midx),1:np));
        end
        for cov_idx=1:n_cov_models
          pidx = (M.patient.cov_idx==cov_idx);
          bagel_Nsil(pidx,c) = sum(M.cov.gene_sil_cov(bagel{c},cov_idx,c),1);
        end
      end
    end

    % bagel flanking
    if P.collapse_categories
      midx = bagel_midx(M.mut.is_flank(bagel_midx));
      if isempty(midx)
        bagel_nflk = zeros(np,1);
      else
        bagel_nflk = as_column(histc(M.mut.pat_idx(midx),1:np));
      end
      bagel_Nflk = zeros(np,1);
      for cov_idx=1:n_cov_models
        pidx = (M.patient.cov_idx==cov_idx);
        bagel_Nflk(pidx) = sum(M.cov.gene_flk_cov_tot(bagel,cov_idx),1);
      end
    elseif ~P.per_category_bagels
      midx = bagel_midx(M.mut.is_flank(bagel_midx));
      bagel_nflk = hist2d_fast(M.mut.pat_idx(midx),M.mut.categ(midx),1,np,1,ncat);
      bagel_Nflk = zeros(np,ncat);
      for cov_idx=1:n_cov_models
        pidx = (M.patient.cov_idx==cov_idx);
        bagel_Nflk(pidx,:) = repmat(shiftdim(sum(M.cov.gene_flk_cov(bagel,cov_idx,:),1))',sum(pidx),1);
      end
    else % P.per_category_bagels
      bagel_nflk = zeros(np,ncat);
      bagel_Nflk = zeros(np,ncat);
      for c=1:ncat
        midx = bagel_midx{c}(M.mut.is_flank(bagel_midx{c}) & M.mut.categ(bagel_midx{c})==c);
        if ~isempty(midx)
          bagel_nflk(:,c) = as_column(histc(M.mut.pat_idx(midx),1:np));
        end
        for cov_idx=1:n_cov_models
          pidx = (M.patient.cov_idx==cov_idx);
          bagel_Nflk(pidx,c) = sum(M.cov.gene_flk_cov(bagel{c},cov_idx,c),1);
        end
      end
    end

  else   % empty bagel
    bagel_nsil=0; bagel_Nsil=0;
    bagel_nflk=0; bagel_Nflk=0;
  end

  % SCALING
  if P.scale_flanking_mutations
    gene_Nflk = gene_Nflk / (globalrate_non/globalrate_flk);
    bagel_Nflk = bagel_Nflk / (globalrate_non/globalrate_flk);  % THIS LINE WAS NOT INCLUDED BEFORE pancan5000/RUN64
    notes_fprintf_once('Scaled flanking mutations (including those of bagel) by adjusting N_flank\n',P.signal);
  end
  if P.scale_silent_mutations
    gene_Nsil = gene_Nsil / (globalrate_non/globalrate_sil);
    bagel_Nsil = bagel_Nsil / (globalrate_non/globalrate_sil);  % THIS LINE WAS NOT INCLUDED BEFORE pancan5000/RUN64
    notes_fprintf_once('Scaled silent mutations (including those of bagel) by adjusting N_flank\n',P.signal);
  end

  % CHOOSE SIGNAL
  if strcmpi(P.signal,'nonsilent')
    n_signal = gene_nnon;
    N_signal = gene_Nnon;
  elseif strcmpi(P.signal,'silent')
    n_signal = gene_nsil;
    N_signal = gene_Nsil;
  elseif strcmpi(P.signal,'flank')||strcmpi(P.signal,'flanking')
    n_signal = gene_nflk;
    N_signal = gene_Nflk;
  else
    error('P.signal = %s not implemented yet',P.signal);
  end

  % CHOOSE TOTAL
  if strcmpi(P.total,'silent+nonsilent+flanking')
    n_total = gene_nnon+gene_nsil+gene_nflk;
    N_total = gene_Nnon+gene_Nsil+gene_Nflk;
  elseif strcmpi(P.total,'silent+nonsilent')
    n_total = gene_nnon+gene_nsil;
    N_total = gene_Nnon+gene_Nsil;
  elseif strcmpi(P.total,'silent+nonsilent+bagel_silent')
    n_total = gene_nnon+gene_nsil+bagel_nsil;
    N_total = gene_Nnon+gene_Nsil+bagel_Nsil;
  elseif strcmpi(P.total,'silent+nonsilent+flanking+bagel_silent+bagel_flanking')
    n_total = gene_nnon+gene_nsil+bagel_nsil+bagel_nflk;
    N_total = gene_Nnon+gene_Nsil+bagel_Nsil+bagel_Nflk;
  else
    error('P.total = %s not implemented yet',P.total);
  end

  if ~isempty(P.special_treatment_for_category)
    ci = P.special_treatment_for_category; if ci<1 || ci>ncat, error('P.special_treatment_for_category out of range'); end
    notes_fprintf_once('Special treatment for category %d:\n',ci);
    notes_fprintf_once('   Signal = gene_S+NS      Total = bagel_S+NS\n');
    n_signal(:,ci) = gene_nnon(:,ci)+gene_nsil(:,ci);
    N_signal(:,ci) = gene_Nnon(:,ci)+gene_Nsil(:,ci);
    n_total(:,ci) = gene_nnon(:,ci)+gene_nsil(:,ci)+bagel_nnon(:,ci)+bagel_nsil(:,ci);
    N_total(:,ci) = gene_Nnon(:,ci)+gene_Nsil(:,ci)+bagel_Nnon(:,ci)+bagel_Nsil(:,ci);
  end

  % make sure values are rounded
  N_total = round(N_total);
  n_total = round(n_total);
  N_signal = round(N_signal);
  n_signal = round(n_signal);

  % make sure we never have N>2*n
  N_total(2*n_total>N_total)=2*n_total(2*n_total>N_total);
  N_signal(2*n_signal>N_signal)=2*n_signal(2*n_signal>N_signal);
  notes_fprintf_once('Making sure we never have N>2*n\n');

  t1=t1+toc(tt); tt = tic;

  % STEP 1
  % for each sample, prioritize mutation categories according to how likely
  % it would be for this gene x sample to have a mutation there by chance.

  P0 = nan(np,ncat);
  for c=1:ncat
    P0(:,c) = my_hygepdf(0,N_total(:,c),n_total(:,c),N_signal(:,c));   % (2D vectorized doesnt work right)
  end
  P0(isnan(P0))=1;
  P1 = 1-P0;

  % determine each patient's priority order of categories (according to P1)
  % left column of "priority" = least extreme category of mutation
  % right column of "priority" = most extreme category of mutation
  [tmp priority] = sort(P1,2,'descend');
  % sort the P arrays to put the columns in least->most priority order
  shft = (priority - repmat(1:ncat,np,1));
  map = reshape(1:(np*ncat),np,ncat);
  newmap = map + shft*np;
  P0 = P0(newmap);
  P1 = P1(newmap);

  t2=t2+toc(tt); tt = tic;

  % STEP 2
  % for each sample, compute probability that it would have been of each degree.
  % where d=0 (no mut) ..... ncat (most extreme mut)
  Pdeg = nan(np,ncat+1);
  for d=0:ncat
    % has to be clear in all categories > d
    Pdeg(:,d+1) = prod(P0(:,d+1:end),2);
    % and nonclear in category d (if d>0)
    if d>0, Pdeg(:,d+1) = Pdeg(:,d+1) .* P1(:,d); end
  end

  %% STEP 2a: calculate score for a sample being of each possible degree
  %% (uses new style, where score = -log10 probability of having a mutation in that category
  %% (zero score for no mutations)
  Sdeg = [zeros(np,1) -log10(P1)];

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
    for d = ncat:-1:1
      c = priority(p,d);
      if n_signal(p,c)>0, degree(p) = d; break; end
    end
    score_obs = score_obs + Sdeg(p,degree(p)+1);
  end

  % for zero score, don't bother doing convolutions
%  if score_obs<=0, G.p(g)=1; continue; end
%  notes_fprintf_once('skipping genes with zero mutations!\n');
  notes_fprintf_once('no longer skipping genes with zero mutations!\n');

  t3=t3+toc(tt); tt = tic;

  % STEP 4
  % compute P value for gene by convolutions
  % --> 95% of the total runtime is in this step (t4/ttot)

  if P.numbins_scheme==1
    numbins = max(10,score_obs*2);
    notes_fprintf_once('binsize fixed at 1; no minimum effect size.  numbins = max(10,score_obs*2)\n');
  elseif P.numbins_scheme==2
    numbins = ceil(score_obs+max(5,0.2*score_obs));
    notes_fprintf_once('binsize fixed at 1; no minimum effect size.  numbins = score_obs+max(5,0.2*score_obs)\n');
  else
    error('unknown P.numbins_scheme');
  end

  % allocate space for convolutions
  H = zeros(numbins,1);
  newH = zeros(numbins,ncols);

  [pmax pmin] = projection_1d_convolutions_fast(Sdeg,Pdeg,score_obs,numbins,H,newH);
  % pmax is the p-value generated by the standard procedure we've always used.
  % pmin is the (better) p-value that goes one discrete step further
  % to solve the problem of discrete statistics, we want to take a randomly chosen position between pmin and pmax.
  G.pmax(g) = pmax;        % for sorting genelist
  G.p(g) = pmin+rand*(pmax-pmin);  % for q-q plot
  G.pmin(g) = pmin;        % for future reference

  t4=t4+toc(tt); tt = tic;
    
  ttot = t1+t2+t3+t4;
    
  progress_fprintf('%5d %12s   %12.10f %12.10f %12.10f   %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n',...
          g,M.gene.name{g},G.pmin(g),G.p(g),G.pmax(g),t1,t2,t3,t4,ttot,ttot/g);

  debug=true;
  if debug
    has_idx=find(any(n_signal,2));
    hasnt_idx = find(~any(n_signal,2));
    pidx = [has_idx(1:min(length(has_idx),13));hasnt_idx(1:min(length(hasnt_idx),2))];
    x=[];x.patient=M.patient.name;x.callscheme=M.patient.callscheme_name;x.priority=priority;
    x.n_signal=n_signal; x.N_signal=round(N_signal);x.n_total=n_total;x.N_total=round(N_total);
    x.degree=degree;x.P0=P0;x.Pdeg=Pdeg;x.Sdeg=Sdeg;
    pr(x,pidx)
%    keyboard
  end

end, fprintf('\n');   % next gene

% FDR
G = sort_struct(G,'pmax');
G.q = calc_fdr_value(G.pmax);

% censor irrelevant columns
if ~P.use_flanking_mutations
  G = rmfield_if_exist(G,{'nflk','Nflk','bagel_nflk','bagel_Nfk'});
end
if ~P.use_bagels
  G = rmfield_if_exist(G,{'bagel_nsil','bagel_nflk','bagel_Nsil','bagel_Nflk'});
end
G = rmfield_if_exist(G,{'bagel'});  % too inconvenient for display

% SAVE+RETURN RESULTS
if exist('outdir','var'), save([outdir '/results.mat'],'G'); end
if exist('outdir','var'), save_struct(G,[outdir '/sig_genes.txt']); end
pr(G,1:min(30,length(gta)))
notes_fprintf('%d genes with q<=0.1\n',sum(G.q<0.105));

if exist('progress_f','var'), fclose(progress_f); end
if exist('notes_f','var'), fclose(notes_f); end




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
