--------------------------------
   SegSeq_running_README.txt
--------------------------------
       TABLE OF CONTENTS
--------------------------------
A. Running SegSeq
B. Calibrating SegSeq parameters
C. Known bugs in SegSeq and suggested improvements
D. Matlab subroutines for SegSeq

[LEGEND]
@@  Bugs
&&  Improvements
**  Notes

-----------------------------------
   A. RUNNING SegSeq
-----------------------------------

(((  A1. BAM FILE PREPROCESSING  )))
** Refer to SegSeq_preprocessing_README.txt for details of
intermediate files

---  A1.1  Syntax for preprocessing script  ---

perl ~/CancerGenomeAnalysis/trunk/matlab/segseq/preprocess/segseq_preprocess_pipeline.pl PARAMETER_LIST

ORDERED PARAMETER_LIST:
   <Normal dir> <Tumor dir> <Sample name>     (Required)
   <Qual cutoff> <Exe dir> <Output dir>       (Optional)

e.g. 
STEP 1
perl
~/CancerGenomeAnalysis/trunk/matlab/segseq/preprocess/segseq_preprocess_pipeline.pl
/xchip/tcga_scratch/ng/OV-0725/wgs/cn/normal_SS/ /xchip/tcga_scratch/ng/OV-0725/wgs/cn/tumor_SS/ test_OV-0725

STEP 2
cd /xchip/tcga_scratch2/ng/test_OV-0725/matfiles
(ls *aligned_reads*mat) -> e.g test_OV-0725_3019J_tumor_normal_aligned_reads_qual20.mat
 
matlab
addpath ~/CancerGenomeAnalysis/trunk/matlab/segseq
SegSeq -s test_output -t test_OV-0725_3019J_tumor_normal_aligned_reads_qual20.mat -wig 10


---  A1.2  Aligned position files  ---
** Replace OV-0751 with the appropriate sample name **

/xchip/tcga_scratch/ng/OV-0751/wgs/cn/normal_SS/chr1.txt .. chr23.txt
/xchip/tcga_scratch/ng/OV-0751/wgs/cn/tumor_SS/chr1.txt .. chr23.txt

Text tab-delimited file with column headers
<LaneID> <ReadID> <Chr> <Start> <End> <Strand> <MappingQuality> <WhichPairmate>


---  A1.3  Lane list files  ---
/xchip/tcga_scratch/ng/OV-0751/wgs/cn/normal_SS/lanelist.txt
/xchip/tcga_scratch/ng/OV-0751/wgs/cn/tumor_SS/lanelist.txt

@@ For each new sample, we need to reconcile the directory tree with
Mike Lawrence's file locations.  I have manually created symbolic
links, e.g. /xchip/tcga_scratch/lawrence/ov/0751/wgs/normal_SS

@@ There is a discrepancy between Mike Lawrence's lanelist.txt file
(0-index) and the LaneID in the aligned position files (1-index)


(((  A2. SEGSEQ SYNTAX  )))

DRIVER FUNCTION:  ~/CancerGenomeAnalysis/trunk/matlab/segseq/SegSeq.m
** Refer to parameter list at the top of the SegSeq.m file for documentation

addpath ~/CancerGenomeAnalysis/trunk/matlab/segseq/
SegSeq -s OV-0751_30W1W -t OV-0751_30W1W_tumor_normal_aligned_paired_reads_qual20.mat
-------------------------------------------------------------------

[For .wig file output with data every (wigStep) reads]
SegSeq -s OV-0751_30W1W -wig 10 -t ...
SegSeq -s OV-0751_30W1W -wig 1  -t ...
-------------------------------------------------------------------

@@ WIG files generated every 1 read (-wig 1) are too big for the
Integrated Genomics Viewer to load or to pre-process


(((  A3. SEGSEQ INPUT  )))

cd /xchip/tcga_scratch2/ng/OV-0751/matfiles/OV-0751

---  A3.1) Matlab file with structures READN, READT  ---
  READN.chr, READN.pos, READN.lane
  READT.chr, READT.pos, READT.lane


---  A3.2)  Optional structure for visualization  ---
  RATIOS.chr, RATIOS.windows, RATIOS.ratios

  % Ratios in 100 kb alignable windows of 36mers can be calculated
  % from this Matlab code snippet:

  addpath ~/CancerGenomeAnalysis/trunk/matlab/segseq/
  load HG18_N36_D2_WINDOWS_100K
  RATIOS=calc_ratios_from_reads( READN, READT, WINDOWS );


(((  A4. SEGSEQ OUTPUT  )))
  SAMPLE_parameters.log      % Log of parameters for sample run
  SAMPLE_parameters.txt      % Text-delimited file with segments
  SAMPLE_parameters.mat      % Matlab structures
   -> SEG                        Final segments after merging
   -> SEG1                       Initial segments before merging
  SAMPLE_parameters.wig      % Local difference statistic in .wig format


(((  A5. VISUALIZING SEGMENTATION RESULTS  )))
  There are three auxiliary plotting scripts for copy-number ratios in
  the ~/CancerGenomeAnalysis/trunk/matlab/segseq

%  Load the current files for the following examples
cd /xchip/tcga_scratch2/ng/OV-0751/matfiles/OV-0751
addpath ~/CancerGenomeAnalysis/trunk/matlab/segseq
load OV-0751_30W1W_tumor_normal_aligned_paired_reads_qual20.mat
load OV-0751_30W1W_noCNV_seg_W_400_initFP_1000_pmerge_1e-20.mat


---  A5.1  Whole genome view  ( plotRatios.R )  ---

plotRatios( RATIOS, SEG );

@@ plotRatios.m only plots chromosomes 1 to 22


---  A5.2  Chromosome view  ( plot_chr_ratios.R )  ---

plot_chr_ratios(RATIOS, 18, 1, 250e6, SEG,'OV-0751_30W1W','g');

% Zoom in on chr18 deletion
axis([32 38 0 1]);


---  A5.3  Cumulative distribution of reads  ( plot_boundaries.m )  ---

plot_boundaries( READN, READT, 18, 34.5, 36, 'OV-0751_30W1W' );



-------------------------------------------------------------------
   B. CALIBRATING SegSeq PARAMETERS
-------------------------------------------------------------------
B1. [SYNTAX]
WARNING: Setting (wigStep) = 1 creates large output files

B2. [DEFAULT PARAMETERS]
-a 1000  Number of genome-wide initial breakpoints from normal replicates
-b 1     Number of genome-wide final segments from normal replicates
-W 400   Number of consecutive normals
-v 1     Filter known copy number variants (use -v 0 to ignore filtering)

[Tuning parameters]
1) p_merge, the genome-wide p-value for the final segmentation, is
determined by the -b 1 parameter on normal replicates.  However, there
is a floor of 1e-20 on p_merge, because G+C bias sometimes makes it
impossible to remove false positive segments in normal replicates.

This p-value cutoff can be made more stringent by invoking the -p parameter

For instance,
  SegSeq -b 1 -s OV-0751 -t OV-0751_tumor_normal_30W1W_aligned_reads.mat
uses a p_merge of 1.05e-20

However, you can set a p_merge override using:
  SegSeq -p 1e-30 -s OV-0751 -t OV-0751_tumor_normal_30W1W_aligned_reads.mat


--------------------------------
   C1. KNOWN BUGS in SegSeq
--------------------------------
1) Hyper-segmentation from lane differences in G+C bias
2) Predicted breakpoints do not agree with funny reads
-> Are duplicate reads and non-unique reads from MAQ alignments
   interfering with local difference statistic?
3) 


--------------------------------------
   C2. OTHER IMPROVEMENTS to SegSeq
--------------------------------------
1) Better choice of lanes for tumor/normal comparisons
  -> Currently just taking all lanes in the flowcell with the median
   number of normal reads
  -> Instead, could choose the 5 tumor lanes and 5 normal lanes with
   the closest G+C content profiles
  





--------------------------------------
   D. MATLAB SUBROUTINES FOR SegSeq.m
--------------------------------------

(PART A) Set segmentation parameters from normal vs. normal comparison
Split reads from half of the normal lanes into N.READN and N.READT
|-> segment_solexa_logratios_normals.m
    ----------------------------------------------------------
     1) CALCULATE BREAKPOINT STATISTICS for normal replicates
    ----------------------------------------------------------
    |->  calc_log_ratio.m                          
         Calculate log difference in copy-number ratios in local windows
         of W reads
    |
    |->  log_normal_approx_for_normratio.m
         Calculate p-values for log difference statistics based on a
         lognormal approximation.  This function uses an fzero solver:
         (diff_log_normaratiopdf_of_exp.m)

    ---------------------------------------------
     2) CALIBRATE (p_bkp) from normal replicates
    ---------------------------------------------
    |-> filter_pval_numFP.m

    -----------------------------------------------
     3) Initialize segments from normal replicates
    -----------------------------------------------
    |-> initialize_seg_logratios.m

    --------------------------------------------------
     4) CALIBRATE p_merge by merging initial segments
        until (numFP) false positives is reached
    -------------------------------------------------
    % NOTE: This is hard-coded to stop at 1e-20 for p-merge %

    while( numBkp > numFinalFP )
    |-> merge_segments_logratio.m
    |-> p_merge = pmerge / 100

(PART B) Repeat segmentation for tumors vs. normals
|-> segment_solexa_logratios.m
    ----------------------------------------------------------
     5) CALCULATE BREAKPOINT STATISTICS for normal replicates
    ----------------------------------------------------------
    |->  calc_log_ratio.m                          
         Calculate log difference in copy-number ratios in local windows
         of W reads
    |
    |->  log_normal_approx_for_normratio.m
         Calculate p-values for log difference statistics based on a
         lognormal approximation.  This function uses an fzero solver:
         (diff_log_normaratiopdf_of_exp.m)

    ---------------------------------------------
     6) Find INITIAL LIST of breakpoints 
    ---------------------------------------------
    |-> filter_pval_v2.m

    -----------------------------------------------
     7) Initialize segments from normal replicates
    -----------------------------------------------
    |-> initialize_seg_logratios.m

    --------------------------------------------------
     8) CALIBRATE p_merge by merging initial segments
        until (numFP) false positives is reached
    -------------------------------------------------
    |-> merge_segments_logratio.m





