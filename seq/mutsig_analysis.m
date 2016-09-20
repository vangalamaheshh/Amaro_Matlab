function M = mutsig_analysis(isetname,outdir,P)

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'build','*required*');  % for COSMIC and clustering analyses
P = impose_default_value(P,'using_scatter_gather',false);
P = impose_default_value(P,'isetname',isetname);

% LOAD DATA
M = load_all_mutation_data(P);

% CREATE OUTPUT DIRECTORY
ensure_dir_exists(outdir);
outstem = [outdir '/' isetname];

% ANALYZE DATA
M = analyze_mutation_data(M,outstem,P);


