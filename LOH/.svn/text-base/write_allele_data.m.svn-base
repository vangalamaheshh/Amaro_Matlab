function write_allele_data(S, bindirname, Lsort, seg_path)

Smin.dat=S.adat(:,:,1);
Smin.chrn=S.chrn;
Smin.chr=S.chr;
Smin.marker=S.marker;
Smin.pos=S.pos;
Smin.sdesc=S.sdesc;
Smax=Smin;
Smax.dat=S.adat(:,:,2);

write_bin([bindirname '_min/'], Smin);
write_bin([bindirname '_max/'], Smax);

write_seg_file([seg_path '.min.ascn.seg.txt'], Smin);
write_seg_file([seg_path '.max.ascn.seg.txt'], Smax);
write_seg_file([seg_path '.loh.seg.txt'], Lsort);
