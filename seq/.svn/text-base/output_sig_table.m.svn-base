function output_sig_table(M,P)

if exist('P','var') && ischar(P)
  tmp = P;
  P = [];
  P.siggene_report_filename = tmp;
end

if ~exist('P','var'), P=[]; end
P=impose_default_value(P,'siggene_report_filename','*required*');

fprintf('Outputting significance table: %s\n',P.siggene_report_filename);

X=[];

if isfield(M.gene,'rank')
  X.rank = M.gene.rank;
else
  X.num = (1:M.ng)';
end

if isfield(M.gene,'genes')   % using genesets
  X.geneset = M.gene.name;
  X.description = M.gene.desc;
  X.genes = cell(M.ng,1); for i=1:M.ng, X.genes{i} = concat(M.gene.genes{i},', '); end
  X.N_genes = as_column(cellfun('length',M.gene.genenos));
  X.mut_tally = M.gene.muttally;
else
  X.gene = M.gene.name;
  if isfield(M.gene,'longname'), X.description = M.gene.longname; end
end

if isfield(M.gene,'N'), X.N = M.gene.N; end
if isfield(M.gene,'n'), X.n = M.gene.n; end
if isfield(M.gene,'npat'), X.npat = M.gene.npat; end
if isfield(M.gene,'nsite'), X.nsite = M.gene.nsite; end
if isfield(M.gene,'nsil'), X.nsil = M.gene.nsil; end

% mutation breakdown by (tumor type) or (category)
flds = grep('^n_',fieldnames(M.gene));
for i=1:length(flds)
  fieldcontents = getfield(M.gene,flds{i});
  tofieldname = flds{i};
  tofieldname = regexprep(tofieldname,'^n_(\D.*)$','$1');
  tofieldname = regexprep(tofieldname,'^n_(\d+)$','n$1');
  X = setfield(X,tofieldname,fieldcontents);
end

% intermediate p_values (if available)
if isfield(M.gene,'pval_classic')
  X.p_classic = ensure_cell(format_number(M.gene.pval_classic,3,8));
  if isfield(M.gene,'pval_classic_lessthan_flag')
    idx = find(M.gene.pval_classic_lessthan_flag);
    X.p_classic(idx) = regexprep(X.p_classic(idx),'^(.*)$','<$1');
  end
end
if isfield(M.gene,'pval_ns_s')
  X.p_ns_s = ensure_cell(format_number(M.gene.pval_ns_s,3,8));
  if isfield(M.gene,'pval_ns_s_lessthan_flag')
    idx = find(M.gene.pval_ns_s_lessthan_flag);
    X.p_ns_s(idx) = regexprep(X.p_ns_s(idx),'^(.*)$','<$1');
  end
end

% MutSig2 p_values (if available)
if isfield(M.gene,'p_clust')
  X.p_clust = ensure_cell(format_number(M.gene.p_clust,3,8));
  if isfield(M.gene,'p_clust_lessthan_flag')
    idx = find(M.gene.p_clust_lessthan_flag);
    X.p_clust(idx) = regexprep(X.p_clust(idx),'^(.*)$','<$1');
  end
end
if isfield(M.gene,'p_cons')
  X.p_cons = ensure_cell(format_number(M.gene.p_cons,3,8));
  if isfield(M.gene,'p_cons_lessthan_flag')
    idx = find(M.gene.p_cons_lessthan_flag);
    X.p_cons(idx) = regexprep(X.p_cons(idx),'^(.*)$','<$1');
  end
end
if isfield(M.gene,'p_joint')
  X.p_joint = ensure_cell(format_number(M.gene.p_joint,3,8));
  if isfield(M.gene,'p_joint_lessthan_flag')
    idx = find(M.gene.p_joint_lessthan_flag);
    X.p_joint(idx) = regexprep(X.p_joint(idx),'^(.*)$','<$1');
  end
end

% final p_value and q_values
X.p = ensure_cell(format_number(M.gene.pval,3,8));
X.q = ensure_cell(format_number(M.gene.qval,3,8));
if isfield(M.gene,'pval_lessthan_flag')
  idx = find(M.gene.pval_lessthan_flag);
  X.p(idx) = regexprep(X.p(idx),'^(.*)$','<$1');
  X.q(idx) = regexprep(X.q(idx),'^(.*)$','<$1');
end

% sort
if isfield(X,'rank')
  X = sort_struct(X,'rank');
end

% save
save_struct(X,P.siggene_report_filename);

fprintf('Done\n');

