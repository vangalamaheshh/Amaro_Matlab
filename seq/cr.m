function cr(R,i)

orf = [];

plusstrand = strcmp(R.strand{i},'+');

if plusstrand, forfrom=1; forstep=+1; forto=R.n_exons(i);
else forfrom=R.n_exons(i); forstep=-1; forto=1;
end

for e=forfrom:forstep:forto
  st = max(R.exon_starts{i}(e),R.code_start(i));
  en = min(R.exon_ends{i}(e),R.code_end(i));
  fprintf('len %d   \tlenmod %d   frame %d   + exon %d  %d-%d\n',...
    length(orf),mod(length(orf),3),R.exon_frames{i}(e),e,st,en);
  if en>=st, orf=[orf genome_region(R.chr(i),st,en)]; end
end



