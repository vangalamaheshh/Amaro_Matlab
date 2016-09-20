function SEG=combine_seg_mat( baseDir, tumorName, suffix, chrs )

SEG.chr = [];
SEG.left = [];
SEG.right = [];
SEG.countN = [];
SEG.countT = [];
SEG.ratios = [];
SEG.diff = [];
SEG.pval = [];


for ci=1:length(chrs)
    c = chrs(ci);

    segfile = [ baseDir tumorName '_chr' num2str(c) '_' suffix ]
    S = load(segfile);

    SEG.chr = [ SEG.chr; S.SEG.chr ];
    SEG.left = [ SEG.left; S.SEG.left ];
    SEG.right = [ SEG.right; S.SEG.right ];
    SEG.countN = [ SEG.countN; S.SEG.countN ];
    SEG.countT = [ SEG.countT; S.SEG.countT ];
    SEG.ratios = [ SEG.ratios; S.SEG.ratios ];
    SEG.diff = [ SEG.diff; S.SEG.diff ];
    SEG.pval = [ SEG.pval; S.SEG.pval ];
end

segfile = [ baseDir '/' tumorName '_' suffix 'seg.txt' ]
write_segments_logratio( segfile, SEG );


matfile = [ baseDir '/' tumorName '_' suffix 'seg.mat' ]
save(matfile, 'SEG');
