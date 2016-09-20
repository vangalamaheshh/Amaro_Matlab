function output_n_and_N_tables(M,P)

if ~exist('P','var'), P=[]; end
P=impose_default_value(P,'filename_stem','table');
P=impose_default_value(P,'sort_table_by_P',true);

% output tables

mutclass = regexprep(M.mutclass,'\+','and');
nc = length(mutclass);
for c=1:nc
  % output n table
  filename = [P.filename_stem '_' mutclass{c} '_mutations.txt'];
  output_table(filename,M,squeeze(M.n_nonsilent(:,c,:)));
  % output N table
  filename = [P.filename_stem '_' mutclass{c} '_coverage.txt'];
  output_table(filename,M,squeeze(M.N_cov(:,c,:)));
end

% output silent n table
  filename = [P.filename_stem '_' mutclass{c} '_silent.txt'];
  output_table(filename,M,squeeze(M.n_silent(:,M.TOT,:)));



function output_table(filename,M,data)

  if P.sort_table_by_P
    Ntot_cov = sum(M.N_cov(:,:,:),3);
    Ntot_tot = Ntot_cov(:,M.TOT);
    ntot_ns = sum(M.n_nonsilent(:,:,:),3);
    ntot_tot = ntot_ns(:,M.TOT);
    [tmp a] = sort(ntot_tot./Ntot_tot, 'descend');
    [tmp b] = sort(M.Prob(a));
    ord = a(b);
  else
    ord = 1:M.ng;
  end

  out = fopen(filename, 'wt');
  fprintf(out, 'gene\tcenter\tp\tq');
  for p=1:M.np
    fprintf(out, '\t');
    if M.patient.hypermutated(p), fprintf(out,'**');
    elseif M.patient.treated(p), fprintf(out,'*');
    end
    fprintf(out, '%s', M.patient.name{p});
  end
  fprintf(out, '\n');
  for i=1:M.ng
    g = ord(i);
    fprintf(out, '%s\t%s\t', M.gene.name{g}, M.gene.center{g});
    if M.Prob(g)==0, fprintf(out, '<1E-11\t<1E-8');
    else fprintf(out, '%s\t%s', format_number(M.Prob(g),3,8), format_number(M.Q(g),3,8));
    end
    for p=1:M.np
      fprintf(out, '\t%d', data(g,p));
    end
    fprintf(out, '\n');
  end
  fclose(out);
end

end
