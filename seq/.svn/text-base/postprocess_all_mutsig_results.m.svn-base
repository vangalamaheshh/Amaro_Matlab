function X = postprocess_all_mutsig_results(P)

if ~exist('P','var'), P=[]; end

P = impose_default_value(P,'runfile','*required*');
P = impose_default_value(P,'mutsigcvdir','*required*');
P = impose_default_value(P,'mutsig2dir','*required*');
P = impose_default_value(P,'maf','*required*'); % only coding mutations needed
P = impose_default_value(P,'outfile','*required*');
P = impose_default_value(P,'cgcfile','/cga/tcga-gsc/home/lawrence/mut/analysis/20130429_pancan/cgc/cancer_gene_census.20130618.ssnv.txt');
P = impose_default_value(P,'p_value_combine_method','fisher');  % 'fisher' or 'truncated_product'

ensure_writeable(P.outfile);

% load runfile
load(P.runfile,'RUN');
X = RUN; clear RUN;
demand_fields(X,{'name','longname','ttype_grep','ttype_grepv','outdir'});

% ensure all necessary files exist
demand_file(P.maf);
for xi=1:slength(X)
  demand_file([P.mutsigcvdir '/' X.outdir{xi} '/results.mat']);
  file = [P.mutsig2dir '/' X.outdir{xi} '/results.mat'];
  if exist(file,'file')
    % OK
  else
    file = [P.mutsig2dir '/' X.name{xi} '/results.mat'];
    if exist(file,'file')
      % OK
    else
      error('mutsig2 file not found for %s',X.name{xi});
    end    
  end
end

% load CGC, KIN, and OG lists, and gene tiers list
C = load_struct(P.cgcfile); %'/cga/tcga-gsc/home/lawrence/mut/analysis/20130429_pancan/cgc/cancer_gene_census.txt');  % 488
C = reorder_struct(C,grepm('M|N|F|S',C.MutationType));               % 162
C.is_dom = grepmi('Dom',C.CancerMolecularGenetics);
C.is_rec = grepmi('Rec',C.CancerMolecularGenetics);
kin = load_kin;
og = load_lines('/cga/tcga-gsc/home/lawrence/mut/analysis/20121031_pancan/overlapping_genes.1.txt');
TIERS=load('/cga/tcga-gsc/home/lawrence/mut/analysis/20130429_pancan/gene_tiers.v1.mat','G'); TIERS=TIERS.G;

% load MAF
fprintf('Loading MAF\n');
M = maftoM(P.maf);
M.mut.is_signal = grepmi('missense|nonsense|splice|ins|del|start|read.?through|non.?stop',M.mut.type);

try

% process each run
fprintf('Processing runs: ');
X.npat = nan(slength(X),1);
X.G = cell(slength(X),1);
for xi=1:slength(X), fprintf('%d/%d ',xi,slength(X));
  cv = load([P.mutsigcvdir '/' X.outdir{xi} '/results.mat'],'G'); cv=cv.G;
  file = [P.mutsig2dir '/' X.outdir{xi} '/results.mat'];
  if exist(file,'file')
    % OK
  else
    file = [P.mutsig2dir '/' X.name{xi} '/results.mat'];
    if exist(file,'file')
      % OK
    else
      error('mutsig2 file not found for %s',X.name{xi});
    end
  end
  m2 = load(file,'G'); m2=m2.G;

  % collect p-values
  if ~isfield(cv,'pmid') && isfield(cv,'p'), cv = rename_field(cv,'p','pmid'); end
  G = rename_field(cv,{'eff','pmin','pmid','pmax'},{'effCV','pminCV','pmidCV','pmaxCV'}); G = rmfield(G,'q');
  G = mapinto(G,m2,'gene','gene',{'nmuts','nperm','effcons','effclust','pcons','pclust','pjoint'});

  % MutSig2:  edit results to weaken instances of nmut=2, pjoint=0, etc.
  %           edit remaining zeros to be 1/nperm
  % MutSigCV: edit zeros to be 1e-16
  pcap = 10.^-(G.nmuts); G.pcons = max(G.pcons,pcap); G.pclust = max(G.pclust,pcap); G.pjoint = max(G.pjoint,pcap);
  pcap = 1./G.nperm;     G.pcons = max(G.pcons,pcap); G.pclust = max(G.pclust,pcap); G.pjoint = max(G.pjoint,pcap);
  pcap = 1e-16;        G.pminCV = max(G.pminCV,pcap); G.pmidCV = max(G.pmidCV,pcap); G.pmaxCV = max(G.pmaxCV,pcap);

  if strcmpi(P.p_value_combine_method,'fisher')
    % Fisher-combine p = pmaxCV and pjoint; calculate q = standard FDR; sort struct
    G.pmid = max(1e-16,fisher_combined_p([G.pmidCV G.pjoint]));
    G.pmax = max(1e-16,fisher_combined_p([G.pmaxCV G.pjoint]));
  elseif strcmpi(P.p_value_combine_method,'truncated_product')
    % Truncated Product Method -combined p = pmaxCV and pjoint
    [tmp G.pmid] = truncated_product_combined_p([G.pmidCV G.pjoint]);
    [G.pmax tmp] = truncated_product_combined_p([G.pmaxCV G.pjoint]);
  else
    error('invalid P.p_value_combine_method');
  end

  G.pmid = max(1e-16,G.pmid);
  G.pmax = max(1e-16,G.pmax);

  % for genes that didn't have MutSig2 run on them, just take the values from CV
  idx = find(isnan(G.pjoint)); G.pmid(idx) = G.pmidCV(idx); G.pmax(idx) = G.pmaxCV(idx);

  % calculate q = standard FDR
  G.q = calc_fdr_value(G.pmax);

  % add npat of this run (includes all mutations)
  % add npat, nsite per gene (includes only coding nonsilent)
  m = M.mut;
  if ~isempty(X.ttype_grep{xi}), m = reorder_struct(m,grepm(X.ttype_grep{xi},m.ttype)); end
  if ~isempty(X.ttype_grepv{xi}), m = reorder_struct_exclude(m,grepm(X.ttype_grepv{xi},m.ttype)); end
  X.npat(xi,1) = length(unique(m.patient));
  m = reorder_struct(m,m.is_signal);
  m.gidx = listmap(m.gene,G.gene);
  if ~isfield(G,'nsite') || ~isfield(G,'npat')
    G.nsite = zeros(slength(G),1);
    G.npat = zeros(slength(G),1);
    for g=1:slength(G)
      idx = find(m.gidx==g); if isempty(idx), continue; end
      G.nsite(g) = length(unique_combos(m.chr(idx),m.start(idx)));
      G.npat(g) = length(unique(m.patient(idx)));
    end
  end
  G = move_field_to_after(G,'npat','nind'); G = move_field_to_after(G,'nsite','npat');

  % anotate with longnames, tiers, is_cgc, is_kin
  % --> blacklist (i.e. set p=q=1) the 'overlapping' problem genes
  G.longname = get_longnames(G.gene); G.tier = mapacross(G.gene,TIERS.gene,TIERS.tier,0);
  G.is_cgc = repmat({'0'},slength(G),1); idx = find(ismember(G.gene,C.Symbol)); G.is_cgc(idx) = repmat({'1'},length(idx),1);
  idx = find(ismember(G.gene,C.Symbol(C.is_dom & ~C.is_rec))); G.is_cgc(idx) = repmat({'1 (Dom)'},length(idx),1);
  idx = find(ismember(G.gene,C.Symbol(~C.is_dom & C.is_rec))); G.is_cgc(idx) = repmat({'1 (Rec)'},length(idx),1);
  G.is_kin = ismember(G.gene,kin);
  % for the overlapping genes, set p=1 and q=1
  G.overlapping = ismember(G.gene,og);
  G = order_fields_first(G,{'gene','longname','tier','is_cgc','is_kin','overlapping'});
  flds = {'pclust','pcons','pjoint','pminCV','pmidCV','pmaxCV','pmid','pmax','q'}; 
  for j=1:length(flds), G.(flds{j})(G.overlapping) = 1; end

  % sort and save
  G = sort_struct(G,{'pmax','npat'},[1 -1]);
  X.G{xi,1} = G;

end,fprintf('\n');   % next run

%%%%%%%%%%%%%%%%%%%%%%%
% UNION
%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Adding UNION\n');
xi = grepvi('union|best_ttype|leaveout|\_|COAD|READ',X.name,1);  % don't count leaveouts or subtypes
U = X.G(xi); for i=1:length(U), U{i}.ttype = repmat({X.name{xi(i)}},slength(U{i}),1); end
U = concat_structs(U); U.q = calc_fdr_value(U.pmax);   % FDR on whole set
% dedup: find the run with the best pmax, and take all the info from that run
U = sort_struct(U,{'pmax','npat'},[1 -1]); [u ui uj] = unique_keepord(U.gene); U = reorder_struct(U,ui);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BEST_TTYPE 
% for red-purple-blue plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Adding BEST_TTYPE\n');
xi = grepvi('union|best_ttype|pan|leaveout|\_|COAD|READ',X.name,1);  % don't count leaveouts or subtypes OR PANCAN
A = X.G(xi); for i=1:length(A), A{i}.ttype = repmat({X.name{xi(i)}},slength(A{i}),1); end
A = concat_structs(A); A.q = calc_fdr_value(A.pmax);   % FDR on whole set
% dedup: find the run with the best pmax, and take all the info from that run
A = sort_struct(A,{'pmax','npat'},[1 -1]); [u ui uj] = unique_keepord(A.gene); A = reorder_struct(A,ui);

%% add union and best_ttype to X
Xu = []; Xu.name = {'union'}; Xu.longname = {'union'};  Xu.npat = max(X.npat);  Xu.G = {U};
Xa = []; Xa.name = {'best_ttype'}; Xa.longname = {'best_ttype'};  Xa.npat = max(X.npat);  Xa.G = {A};
X = concat_structs_keep_all_fields({Xu,Xa,X});

%X = order_fields_first(X,{'name','npat','longname','ttype_grep','ttype_grepv','outdir','G','ncons','nclust','njoint2','nCV','njoint3'});
X = order_fields_first(X,{'name','npat','longname','ttype_grep','ttype_grepv','outdir','G'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COUNT NUMBERS OF SIGNIFICANT GENES IN EACH RUN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X.ncons = nan(slength(X),1);
X.nclust = nan(slength(X),1);
X.njoint2 = nan(slength(X),1);
X.nCV = nan(slength(X),1);
X.njoint3 = nan(slength(X),1);

for xi=1:slength(X)
  if xi>2    
    q = calc_fdr_value(X.G{xi}.pcons); X.ncons(xi)=sum(q<=0.1);
    q = calc_fdr_value(X.G{xi}.pclust); X.nclust(xi)=sum(q<=0.1);
    q = calc_fdr_value(X.G{xi}.pjoint); X.njoint2(xi)=sum(q<=0.1);
    q = calc_fdr_value(X.G{xi}.pmaxCV); X.nCV(xi)=sum(q<=0.1);
  end
  X.njoint3(xi,1)=sum(X.G{xi}.q<=0.1);
end

% PRINT SUMMARY OF GENECOUNTS
pr(X,{'name','npat','ncons','nclust','njoint2','nCV','njoint3'});

% SAVE OUTPUT
save(P.outfile,'X');

catch me, excuse(me); end


