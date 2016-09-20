function output_gene_sample_coverage_table(M,outname,P)

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'simplified_gene_sample_coverage_table',true);

require_fields(M,{'gene','patient','np','N_terr','N_cov','TOT'});
require_fields(M.patient,{'name'});

fprintf('Outputting gene_sample_coverage_table: %s\n',outname)

if P.simplified_gene_sample_coverage_table

  X=M.gene;
  fieldname = genfieldname(regexprep(M.patient.name,'-','_'));
  for p=1:M.np
    X = setfield(X,fieldname{p},M.N_cov(:,M.TOT,p));
  end
  save_struct(X,outname);

else

  require_fields(M,{'ng'});
  require_fields(M.gene,{'name','chr','start','end'});
  M.gene.chr = convert_chr(M.gene.chr);
  M.gene = make_numeric(M.gene,{'start','end','len'});

  out = fopen(outname,'wt');
  fprintf(out,'gene\tchr\tmin\tmax\tlen');
  for p=1:M.np, fprintf(out,'\t%s',M.patient.name{p}); end
  fprintf(out,'\n');
  for g=1:M.ng, if ~mod(g,1000), fprintf('%d/%d ',g,M.ng); end
    fprintf(out,'%s\t%d\t%d\t%d\t',M.gene.name{g},M.gene.chr(g),M.gene.start(g),M.gene.end(g));
    for c=1:M.TOT
      fprintf(out,'%d',M.N_terr(g,c));
      if c<M.TOT, fprintf(out,','); end
    end
    for p=1:M.np
      fprintf(out,'\t');
      for c=1:M.TOT
        fprintf(out,'%d',M.N_cov(g,c,p));
        if c<M.TOT, fprintf(out,','); end
      end
    end
    fprintf(out,'\n');
  end, fprintf('\n');
  fclose(out);

end
