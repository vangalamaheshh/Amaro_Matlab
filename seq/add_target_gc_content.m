function add_target_gc_content(infile,outfile,build)

if ~exist('build','var')
  build = 'hg18';
  fprintf('Assuming build hg18\n');
end

try

%X = load_struct_noheader(infile);
X = rename_fields(load_struct_noheader(infile),...
  {'col1','col2','col3','col4'},{'gene','chr','start','end'});
X = make_all_possible_numeric(X);
X.chr = convert_chr(X.chr);
nx = slength(X);
X.gc = nan(nx,1);
for chr=1:24
  fprintf('Chr %d\n',chr);
  ref = upper(genome_region(chr,1,inf,build));
  len = length(ref);
  gc = nan(len,1);
  gc(ref=='A' | ref=='T') = 0;
  gc(ref=='C' | ref=='G') = 1;
  idx = find(X.chr==chr);
  for j=1:length(idx), i=idx(j);
    st = max(1,X.start(i));
    en = min(len,X.end(i));
    X.gc(i) = nanmean(gc(st:en));
  end
end


save_struct(X,outfile,'no_headers');

catch me, excuse(me); end
