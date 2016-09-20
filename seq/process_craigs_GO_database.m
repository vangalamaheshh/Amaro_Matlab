function process_craigs_GO_database(infile,outfile)

fprintf('Loading input file\n');
GO = load(infile);
ngo = length(GO.gene_lists);

fprintf('Writing output file\n');
f = fopen(outfile,'wt');
for i=1:ngo, if ~mod(i,1000), fprintf('%d/%d ',i,ngo); end
  fprintf(f,'%s\t%s',GO.Ds{i,1},GO.Ds{i,2});
  for j=1:length(GO.gene_lists{i})
    fprintf(f,'\t%s',GO.gene_lists{i}{j});
  end
  fprintf(f,'\n');
end, fprintf('\n');
fclose(f);


