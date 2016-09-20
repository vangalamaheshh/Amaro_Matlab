#! /bin/csh -f
# code from Gordon Saksena 2009-08-03

#do preprocess and segseq on all lanes on one sample

setenv SAMP $1

cd ~gsaksena/CancerGenomeAnalysis/trunk/segseq/preprocess3

matlab -nodisplay -r "segseq_preprocess3('/xchip/tcga_scratch/ng/$SAMP/wgs/cn/normal_SS', '/xchip/tcga_scratch/ng/$SAMP/wgs/cn/tumor_SS', '$SAMP', 0);exit;"
# used to be matlab -nodesktop

cd /xchip/tcga_scratch2/ng/$SAMP/matfiles

foreach FLOWCELLFILE (*tumor_normal_aligned_paired_reads*mat)

setenv FLOWCELL `echo $FLOWCELLFILE | sed 's/_tumor.*$//g' | sed 's/_medianfc//g'| sed 's/_nonmedfc//g' `

matlab -nodisplay -r "addpath ~gsaksena/CancerGenomeAnalysis/trunk/segseq; try SegSeq -s $FLOWCELL -t $FLOWCELLFILE -wig 10);catch, sin(1); end; exit;"
# used to be matlab -nodesktop

end

