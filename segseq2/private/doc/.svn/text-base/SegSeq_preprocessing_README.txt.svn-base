--------------------------------
       TABLE OF CONTENTS
--------------------------------
A. File format & directory summary
B. Known bugs in SegSeq and suggested improvements
C. Pre-processing pipeline


----------------------------------------
   A. FILE FORMAT & DIRECTORY SUMMARY
----------------------------------------

1) Need to create symbolic links to Mike Lawrence's directories:
/xchip/tcga_scratch/ng/OV-0725/wgs/cn/normal_SS  ->
/xchip/tcga_scratch/lawrence/ov/0725/wgs/normal_SS

/xchip/tcga_scratch/ng/OV-0725/wgs/cn/tumor_SS -> 
/xchip/tcga_scratch/lawrence/ov/0725/wgs/tumor_SS

2) Temporary Matlab structures
/xchip/tcga_scratch2/ng/OV-0725/dat/

3) Aligned reads input for SegSeq
/xchip/tcga_scratch2/ng/OV-0725/matfiles/OV-0725_tumor_normal_aligned_paired_reads_qual20.mat

4) SegSeq output
[Segment file - TXT]
/xchip/tcga_scratch2/ng/OV-0725/matfiles/OV-0725_302TV_seg_W_400_initFP_1000_pmerge_4.2128e-32.txt

[Segment file - MAT]
/xchip/tcga_scratch2/ng/OV-0725/matfiles/OV-0725_302TV_seg_W_400_initFP_1000_pmerge_4.2128e-32.mat

[Local difference statistic - WIG]
/xchip/tcga_scratch2/ng/OV-0725/matfiles/OV-0725_302TV_seg_W_400_initFP_1000.wig



--------------------------------
   B. KNOWN BUGS in SegSeq
--------------------------------
1) Copy number variants not being filtered
2) p_merge is hard-coded to 1e-20; does not converge with some tumor
   -> Copy number variants?
3) Hyper-segmentation from lane differences in G+C bias
4) Predicted breakpoints do not agree with funny reads
5) Are duplicate reads and non-unique reads interfering with local
   difference statistic?
6) 


--------------------------------
   B. SUGGESTED IMPROVEMENTS to pre-processing pipeline
--------------------------------
1) Automate symbolic links to Mike Lawrence's directories
2) Delete temporary Matlab structures:
   /xchip/tcga_scratch2/ng/[SAMPLE]/dat/chr*mat
   /xchip/tcga_scratch2/ng/[SAMPLE]/matfiles/[SAMPLE]_normal_aligned_paired_reads_qual20.mat
   /xchip/tcga_scratch2/ng/[SAMPLE]/matfiles/[SAMPLE]_tumor_aligned_paired_reads_qual20.mat



--------------------------------
   C. PRE-PROCESSING PIPELINE
--------------------------------
PRE-PROCESSING STEPS
1) Obtain aligned reads from BAM files
2) Convert Mike Lawrence's text files to sorted Matlab structures
3) Collect reads across chromosomes
4) Correct G+C content
5) Split reads by flowcells; only save flowcell with the median number
   of reads
--------------------------------

perl ~/CancerGenomeAnalysis/sandbox/dchiang/segseq_preprocess_pipeline.pl
   <Normal dir> <Tumor dir> <Sample name> [ <Qual cutoff> <Exe dir>
   <Output dir> ]

Optional paramters in [ ]

e.g. perl
~/CancerGenomeAnalysis/sandbox/dchiang/segseq_preprocess_pipeline.pl
/xchip/tcga_scratch/ng/OV-0751/normal_SS/
/xchip/tcga_scratch/ng/OV-0751/tumor_SS/ OV-0751



C.1.1) ALIGNED READ POSITIONS
--------------------------
[Created by]   /xchip/tcga/gbm/analysis/lawrence/SegSeqPreprocess.java

/xchip/tcga_scratch/lawrence/[SAMPLE_NAME]/normal_SS/chr1.txt ... chr24.txt
/xchip/tcga_scratch/lawrence/[SAMPLE_NAME]/tumor_SS/chr1.txt ... chr24.txt

Text tab-delimited file with column headers
<LaneID> <ReadID> <Chr> <Start> <End> <Strand> <MappingQuality> <WhichPairmate>


C.1.2) LANE/FLOWCELL CORRESPONDENCE
------------------------------------
[Created by]    /xchip/tcga/gbm/analysis/lawrence/SegSeqMakeLaneList.java

/xchip/tcga_scratch/lawrence/[SAMPLE_NAME]/normal_SS/lanelist.txt
/xchip/tcga_scratch/lawrence/[SAMPLE_NAME]/tumor_SS/lanelist.txt

Text tab-delimited file with column headers
<LaneID> <Read group>

*** KNOWN ISSUE ***  
The LaneID do not match between file types:
<LaneID> in the lanelist.txt files uses 0-indexing
<LaneID> in the chr#.txt files uses 1-indexing



C.2.1) BINARY MATLAB STRUCTURES FOR ALIGNED READS
----------------------------------------------------
[Created by]    convert_aligned_pairs_to_mat
                -> Compiled Matlab code


