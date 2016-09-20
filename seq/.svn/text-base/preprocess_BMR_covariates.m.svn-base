function [cov_matrix_tot_nopats_nocats, mut_matrix_tot_nopats_nocats, E1, M2] =  preprocess_BMR_covariates(M, BMR_covariates_file)

%E = load_struct(BMR_covariates_file, ['%s', '%f' ]);
%M = load_mutsig(mutsig_outdir, 'LUSC', 'hg19');

E = load_struct_noheader(BMR_covariates_file);
if strcmpi(E.col1{1},'gene') || strcmpi(E.col1{1},'name')
  E = reorder_struct(E,2:slength(E));
end
E = make_numeric(E,'col2');
E = rename_fields(E,{'col1','col2'},{'gene','expr'});

%% TO DO->>>>>remove dependence on "median" field in load_all_mutation_data
%%      ->>>>>make compatible with multiple columns

%% Intersect gene names 

[inters midx] = intersect(M.gene.name, E.gene);
[inters eidx] = intersect(E.gene, M.gene.name);
E = reorder_struct(E, eidx);
M2 = reorder_struct(M.gene, midx);

%% Reorganize

mut_matrix_tot = M.n_silent + M.n_nonsilent;
mut_matrix_ns = M.n_nonsilent;
mut_matrix_s = M.n_silent ;
mut_matrix_tot = mut_matrix_tot(midx, :, :);
mut_matrix_ns = mut_matrix_ns(midx, :, :);
mut_matrix_s = mut_matrix_s(midx, :, :);

cov_matrix_tot = M.N_sil_cov + M.N_non_cov;
cov_matrix_s = M.N_sil_cov;
cov_matrix_ns = M.N_non_cov;
cov_matrix_tot = cov_matrix_tot(midx, :, :);
cov_matrix_s = cov_matrix_s(midx, :, :);
cov_matrix_ns = cov_matrix_ns(midx, :, :);

% sort by name
[e ei] = sort(E.gene);
[m2 m2i] = sort(M2.name);
E1 = reorder_struct(E, ei);
M2 = reorder_struct(M2, m2i);
mut_matrix_tot = mut_matrix_tot(m2i, :, :);
mut_matrix_s = mut_matrix_s(m2i, :, :);
mut_matrix_ns = mut_matrix_ns(m2i, :, :);
cov_matrix_tot = cov_matrix_tot(m2i, :, :);
cov_matrix_s = cov_matrix_s(m2i, :, :);
cov_matrix_ns = cov_matrix_ns(m2i, :, :);


% sort by expression
[e2 ei2] = sort(E.expr);
E1 = reorder_struct(E1, ei2);
M2 = reorder_struct(M2, ei2);
mut_matrix_tot = mut_matrix_tot(ei2, :, :);
mut_matrix_s = mut_matrix_s(ei2, :, :);
mut_matrix_ns = mut_matrix_ns(ei2, :, :);
cov_matrix_tot = cov_matrix_tot(ei2, :, :);
cov_matrix_s = cov_matrix_s(ei2, :, :);
cov_matrix_ns = cov_matrix_ns(ei2, :, :);

% collapse matrices
mut_matrix_tot_nopats = sum(mut_matrix_tot, 3);
mut_matrix_s_nopats = sum(mut_matrix_s, 3);
mut_matrix_ns_nopats = sum(mut_matrix_ns, 3);
cov_matrix_tot_nopats = sum(cov_matrix_tot, 3);
cov_matrix_s_nopats = sum(cov_matrix_s, 3);
cov_matrix_ns_nopats = sum(cov_matrix_ns, 3);

mut_matrix_tot_nopats_nocats = mut_matrix_tot_nopats(:,7);
mut_matrix_s_nopats_nocats = mut_matrix_s_nopats(:,7);
mut_matrix_ns_nopats_nocats = mut_matrix_ns_nopats(:,7);
cov_matrix_tot_nopats_nocats = cov_matrix_tot_nopats(:,7);
cov_matrix_s_nopats_nocats = cov_matrix_s_nopats(:,7);
cov_matrix_ns_nopats_nocats = cov_matrix_ns_nopats(:,7);

mut_matrix_tot_nocats = mut_matrix_tot(:,7);
mut_matrix_s_nocats = mut_matrix_s(:,7);
mut_matrix_ns_nocats = mut_matrix_ns(:,7);
cov_matrix_tot_nocats = cov_matrix_tot(:,7);
cov_matrix_s_nocats = cov_matrix_s(:,7);
cov_matrix_ns_nocats = cov_matrix_ns(:,7);

