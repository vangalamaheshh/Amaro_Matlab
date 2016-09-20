function TSP = load_TSP
%
% load_TSP(P)
%

TSP = load_mutdata2('TSP',188,623);

P=[];
P.use_partitioned_rollup_coverage = true;
P.coverage_replace_cutoff = 0;
TSP = choose_mutations(TSP,P);
TSP = build_n_and_N_tables(TSP,P);
TSP = replace_missing_coverage(TSP,P);

