function G = calcsig(G,P)

demand_field(G,{'name','n_work','N_work'});

P.gene_names = G.name;

G.ntot = sum(G.n_work,3);
G.Ntot = sum(G.N_work(:,end,:),3);

if isfield(G,'X_work') && isfield(G,'x_work') && ~isfield(G,'mu_work')
  G.p = calculate_significance(G.N_work,G.n_work,{G.X_work,G.x_work},P);
elseif ~isfield(G,'X_work') && ~isfield(G,'x_work') && isfield(G,'mu_work')
  G.p = calculate_significance(G.N_work,G.n_work,G.mu_work,P);
else
  error('G needs either mu_work OR X_work, x_work');
end

G.q = calc_fdr_value(G.p);

G = sort_struct(G,'p');

