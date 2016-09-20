function refs = find_references_that_overlap_target(D,targ)
ridx = find(strcmp(D.chrom,targ.chr) & D.txStart <= max(targ.region_ends) & D.txEnd >= min(targ.region_starts));
refs=cell(length(ridx),1);
for i=1:length(ridx)
  refs{i}.build = targ.build;
  refs{i}.chr = D.chrom{ridx(i)};
  refs{i}.strand = D.strand{ridx(i)};
  refs{i}.coding_start = D.cdsStart(ridx(i));
  refs{i}.coding_end = D.cdsEnd(ridx(i));
  refs{i}.exon_starts = D.exonStarts{ridx(i)};
  refs{i}.exon_ends = D.exonEnds{ridx(i)};
end
