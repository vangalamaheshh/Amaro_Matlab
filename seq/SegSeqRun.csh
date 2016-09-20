#! /bin/csh -f
# code from Gordon Saksena 2009-08-03

#do preprocess and segseq on all lanes on one sample

setenv NORM_SS $1
setenv TUM_SS $2
setenv SAMPLE_NAME $3
setenv OFFSET_MISMATCH $4
setenv QUAL_CUTOFF $5
setenv EXE_DIR $6
setenv BASE_DIR $7
setenv MAX_CHR $8
setenv OUTPUT_ALL_FLOWCELLS $9

cd ~gsaksena/CancerGenomeAnalysis/trunk/segseq/preprocess3

reuse matlab
matlab -nodisplay -r "segseq_preprocess3('$1', '$2', '$3', '$4', '$5', '$6', '$7', '$8', '$9');exit;"

cd $BASE_DIR/matfiles

foreach FLOWCELLFILE (*tumor_normal_aligned_paired_reads*mat)

setenv FLOWCELL `echo $FLOWCELLFILE | sed 's/_tumor.*$//g' | sed 's/_medianfc//g'| sed 's/_nonmedfc//g' `

matlab -nodisplay -r "addpath ~gsaksena/CancerGenomeAnalysis/trunk/segseq; try SegSeq -s $FLOWCELL -t $FLOWCELLFILE -wig 10);catch, sin(1); end; exit;"

end

