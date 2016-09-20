function Msil = prepare_for_silent_analysis(M,P)

if ~exist('P','var'), P=[]; end

Msil = M;
Msil.use = M.use_silent;
Msil.use_nonsilent = M.use_silent;
Msil = rmfield(Msil,'use_silent');

if isfield(M,'N_sil_cov')
  Msil.N_non_cov = M.N_sil_cov;
  Msil = rmfield(Msil,'N_sil_cov');
end

Msil.n_nonsilent = M.n_silent;

if isfield(M,'n_nonsilent_ignoring_null_categ')
  Msil.n_nonsilent_ignoring_null_categ = M.n_silent;
end

Msil = rmfield(Msil,'n_silent');
Msil.note1 = 'FOR SILENT P-VALUE CALC';



