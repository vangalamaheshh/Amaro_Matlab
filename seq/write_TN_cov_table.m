function write_TN_cov_table(C,outstem)

try

for i=1:2
  if i==1, tn='tumor'; else tn='normal'; end
  outfile = [outstem tn '_cov.txt'];
  fprintf('Writing %s    ... ', outfile);
  out = fopen(outfile,'wt');
  fprintf(out,'gene\tchr\tstart\tend\tlen\tgc');
  for s=1:C.ns,fprintf(out,'\t%s',C.sample.med{s});end
  fprintf(out,'\n');
  for g=1:C.ng
    fprintf(out,'%s\t%d\t%d\t%d\t%d\t%0.2f',C.gene.name{g},C.gene.chr(g),...
       C.gene.start(g),C.gene.end(g),C.gene.len(g),C.gene.gc(g));
    for s=1:C.ns
      if i==1, ct=C.gtumcov(g,s); else ct=C.gnormcov(g,s); end
      fprintf(out,'\t%d',ct);
    end
    fprintf(out,'\n');
  end
  fclose(out);
  fprintf('done.\n');
end

catch me; excuse(me); end
