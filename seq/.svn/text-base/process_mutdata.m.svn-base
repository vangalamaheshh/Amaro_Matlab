function M = process_mutdata(M,P)

M = choose_mutations(M,P);
M = build_n_and_N_tables(M,P);
M = process_3N_breakdown_stats(M);
M = replace_missing_coverage(M,P);
M = define_indel_bases_at_risk(M);
M = remove_impossible_mutations(M);
M = collapse_mutation_categories(M,P);
M = analyze_mutation_rates(M,P);

