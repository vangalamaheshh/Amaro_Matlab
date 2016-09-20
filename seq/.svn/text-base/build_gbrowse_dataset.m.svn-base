function build_gbrowse_dataset(P)

if ~exist('P','var'), P=[]; end

P = impose_default_value(P,'runfile','*required*');
P = impose_default_value(P,'datasetdir','*required*');
P = impose_default_value(P,'title','*required*');
P = impose_default_value(P,'maf','*required*'); % should be Oncotator MAF for most usefulness to website users
P = impose_default_value(P,'cnfile',''); % for GISTIC results
P = impose_default_value(P,'cgcfile','/cga/tcga-gsc/home/lawrence/mut/analysis/20130429_pancan/cgc/cancer_gene_census.20130618.ssnv.txt');

fprintf('Loading runfile.\n');
load(P.runfile,'X');

demand_file(P.maf);
if ~isempty(P.cnfile), demand_file(P.cnfile); end
demand_file(P.cgcfile);
ede(P.datasetdir);

% save ttypes list
fprintf('Writing runs/tumortypes list.\n');
ttypes = keep_fields(X,{'name','longname','npat'});
ttypes.num = (1:slength(ttypes))';
ttypes = order_fields_first(ttypes,'num');
save_struct(ttypes,[P.datasetdir '/ttypes.txt']);

% save title
save_textfile(P.title,[P.datasetdir '/title.txt']);

% for each tumotype/run, make siggene table, sorted by min of the three q-values
fprintf('Writing per-tumortype genelists... ');
outdir = [P.datasetdir '/siggenes']; ede(outdir);
for i=1:slength(X),fprintf('%d/%d ',i,slength(X));
  x = X.G{i};
  x = keep_fields_if_exist(x,{'gene','longname','tier','is_cgc','is_kin','codelen','exprmax','rt','hiC','paz','nstrikes','nnei',...
                      'nncd','nsyn','nmis','nnon','nspl','nind','npat','nsite','pmaxCV','pcons','pclust','pmax','q','ttype'});
  x = rename_fields(x,{'pmaxCV','pmax'},{'pCV','p'});
  x = sort_struct(x,{'p','npat'},[1 -1]);
  save_struct(x,[outdir '/ttype' num2str(i) '.sig_genes.txt']);
end,fprintf('\n');

% make CGC list
fprintf('Writing CGC list.\n');
C = load_struct(P.cgcfile);
C = reorder_struct(C,grepm('M|N|F|S',C.MutationType));  % (this step may have already been done)
C.domres = regexprep(C.CancerMolecularGenetics,'^(.){1,8}.*$','$1');
C.domres = regexprep(C.domres,'\"','');
C.ttype = regexprep(C.TumourTypesSomaticMutations,'^(.){1,40}.*$','$1');
C.ttype = regexprep(C.ttype,'\"','');
if isfield(C,'is2004') && ~isfield(C,'in2004'), C = rename_field(C,'is2004','in2004'); end
if isfield(C,'in_futreal2004') && ~isfield(C,'in2004'), C = rename_field(C,'in_futreal2004','in2004'); end
C = keep_fields_if_exist(C,{'domres','ttype','in2004','gene','ttype_included'});
if isfield(C,'ttype_included'), C = move_field_to_after(C,'ttype_included','ttype'); C = rename_field(C,'ttype_included','ttype_incl'); end
if isfield(C,'ttype'), C = rename_field(C,'ttype','cgc_ttype'); end
tti = find(strcmp(X.name,'union'));
C = mapinto(C,X.G{tti},'gene','gene',{'codelen','nncd','nsyn','nmis','nnon','nspl','nind','npat','nsite','pmaxCV','pcons','pclust','pmax','q','ttype'},...
            {'codelen','nncd','nsyn','nmis','nnon','nspl','nind','npat','nsite','pCV','pcons','pclust','p','q','ttype'});
C = sort_struct(C,{'p','npat'},[1 -1]);
save_struct(C,[P.datasetdir '/cgc_results.txt']);
save([P.datasetdir '/cgc_results.mat'],'C');
% also save version restricted to tumor types we have
if isfield(C,'ttype_incl')
  C1 = reorder_struct(C,strcmp(C.ttype_incl,'1'));
  save_struct(C1,[P.datasetdir '/cgc_results.ttype_included.txt']);
end

% for each gene, make per-tumortype table, and info file
fprintf('Writing per-gene per-tumorttype tables... ');
outdir = [P.datasetdir '/genetables']; ede(outdir);
f1 = {'gene','longname','tier','is_cgc','is_kin','codelen','exprmax','rt','hiC','paz','nstrikes'};
G = keep_fields_if_exist(X.G{1},f1);
for i=1:slength(G), if ~mod(i,100), fprintf('%d/%d ',i,slength(G)); end
  save_struct(reorder_struct(G,i),[outdir '/' G.gene{i} '.info.txt']);
  TT=[]; TT.ttype = X.longname; TT.ttype_npat = X.npat;
  f = {'nncd','nsyn','nmis','nnon','nspl','nind','npat','nsite','pmaxCV','pcons','pclust','pmax','q'};
  gidx = nan(slength(TT),1);
  for tti=1:slength(X)
    gidx(tti) = find(strcmp(G.gene{i},X.G{tti}.gene));
  end
  for fi=1:length(f)
    TT.(f{fi}) = nan(slength(TT),1);
    for tti=1:slength(X)
      TT.(f{fi})(tti) = nansub(X.G{tti}.(f{fi}),gidx(tti));
    end
  end
  TT = rename_field(TT,{'pmaxCV','pmax'},{'pCV','p'});
  save_struct(TT,[outdir '/' G.gene{i} '.by_ttype.txt']);
end, fprintf('\n');

% for each gene, make gistic results file
load(P.cnfile,'cn');
G = reorder_struct(cn.gene,ismember(cn.gene.name,X.G{1}.gene));
fprintf('Writing GISTIC results files... ');
outdir = [P.datasetdir '/gistic']; ede(outdir);
for i=1:slength(G), if ~mod(i,100), fprintf('%d/%d ',i,slength(G)); end
  tmp=[];
  tmp.amps = G.gistic_amps(i);
  tmp.dels = G.gistic_dels(i);
  save_struct(tmp,[outdir '/' G.name{i} '.gistic.txt']);
end, fprintf('\n');

% make per-gene MAF files
fprintf('Loading MAF.\n');
m = load_struct(P.maf);
fprintf('Writing per-gene MAF files... ');
outdir = [P.datasetdir '/genemafs']; ede(outdir);
[u ui uj] = unique(m.Hugo_Symbol);
for i=1:length(u), if ~mod(i,100), fprintf('%d/%d ',i,length(u)); end
  if strcmp(u{i},'')||strcmpi(u{i},'Unknown')||strcmp(u{i},'?')||strcmpi(u{i},'IGR')||strcmp(u{i},'-'), continue; end
  save_struct(reorder_struct(m,uj==i),[outdir '/' u{i} '.maf']);
end, fprintf('\n');

