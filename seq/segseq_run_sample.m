function segseq_run_sample(normal_ss,tumor_ss,sample_name,offset_mismatch,qual_cutoff,exe_dir,base_dir,max_chr, ...
                           output_all_flowcells)

addpath /xchip/cga2/lawrence/cga/trunk/segseq
addpath /xchip/cga2/lawrence/cga/trunk/segseq/preprocess3
addpath /xchip/cga2/lawrence/cga/trunk/segseq/SegSeqSplitFlowCells
addpath /xchip/cga2/lawrence/cga/trunk/segseq/utilities
addpath /xchip/cga2/lawrence/cga/trunk/segseq/SegSeqRun
addpath /xchip/cga2/lawrence/cga/trunk/segseq/qltout

segseq_preprocess3(normal_ss,tumor_ss,sample_name,offset_mismatch,qual_cutoff,exe_dir,base_dir,max_chr, ...
                           output_all_flowcells)

segseq_run_each_flowcell(base_dir);
