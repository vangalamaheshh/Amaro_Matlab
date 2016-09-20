To prepare flat files for the website run make_full_from_segfiles 
(which prepares and executes all the necessary GISTIC runs), then run 
TCGA_website_flat_file (runs the analysis). The output path has to be changed 
in both modules.

The definitions of cancer types and their hierarchical relation is defined in
tcga_cancer_types. The location of input segfiles is in tcgascape_segfiles 
(currently hardwired to genepattern outputs from a Jan 14 2011 run).

The other modules are more or less generic:

count_cancer_types.m - counts number of samples for each cancer type
run_gistic_on_all.m - runs GISTIC on every D.<type>.mat in a given directory
load_D_tree.m - load all the segmented input ino one structure
clean_gistic_input.m - CNV, NaN, non-autosomal chromosome removal + smoothing
splitsave_D.m - save a D struct off as hierarchical types.
index_cancer_types.m - sample indices for all cancer types 
leaf_cancer_types.m - indicator if type is independent (leaf) type or not

merge_Ds.m, make_IGV_inputs.m, unify_tcgascape_Ds.m, make_master_D.m, and 
unipos_D.m are no longer used since our inputs are now segfiles and not GISTIC 
outputs.

Steve Schumacher 2011.01.31
