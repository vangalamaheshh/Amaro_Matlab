function M = calcsig_v1(M,P)

demand_field(M.gene,{'name','n_work','N_work'});

P.gene_names = M.gene.name;

M.gene.ntot = sum(M.gene.n_work,3);
M.gene.Ntot = sum(M.gene.N_work(:,end,:),3);

M.gene = rmfield_if_exist(M.gene,{'p','q','rank'});
M.gene = orderfields_first(M.gene,{'name'});

if isfield(M.gene,'X_work') && isfield(M.gene,'x_work') && ~isfield(M.gene,'mu_work')
  M.gene.p = calculate_significance(M.gene.N_work,M.gene.n_work,{M.gene.X_work,M.gene.x_work},P);
elseif ~isfield(M.gene,'X_work') && ~isfield(M.gene,'x_work') && isfield(M.gene,'mu_work')
  M.gene.p = calculate_significance(M.gene.N_work,M.gene.n_work,M.gene.mu_work,P);
else
  error('M.gene needs either mu_work OR X_work, x_work');
end

M.gene.q = calc_fdr_value(M.gene.p);

if isfield(M.gene,'effect')
  [tmp ord] = sort_struct(M.gene,{'p','effect'},[1 -1]);
else
  [tmp ord] = sort_struct(M.gene,'p');
end

[tmp M.gene.rank] = sort(ord);



