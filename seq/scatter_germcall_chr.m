function scatter_germcall_chr(boomdir,outdir,chr)

if chr<1 || chr>24, error('chr should be 1-24'); end

chrlen=load_chrlen;
chrlen=chrlen(chr);
params = [];
params.lod_type = 'other_vs_ref';
params.lod_cutoff = 5;
params.lod_cutoff_2 = 0;
params.coverage_cutoff = 15;

M={};idx=1;
for first=1:1e6:chrlen
  last=min(chrlen,first+1e6-1);
  fprintf('\nchr%d:%d-%d\n',chr,first,last);
  try
    M{idx} = germcall(boomdir,chr,first,last,params);
  catch me
    M{idx} = zeros(0,9);
  end
  idx=idx+1;
end
M = cat(1,M{:});

if ~exist(outdir,'dir'), mkdir(outdir); end
T = print_mutation_list(M);
save_textfile(T,[outdir '/chr' num2str(chr) '.txt']);
