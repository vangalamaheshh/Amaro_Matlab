function scatter_somcall_chr(t_boomdir,n_boomdir,outdir,chr,params)
% scatter_somcall_chr(t_boomdir,n_boomdir,outdir,chr,params)
%
% helper function to scatter_somcall
% for calling somatic mutations in exome
%
% Mike Lawrence 2009

if ~exist('params','var'), params=[]; end
if ischar(params)
  tmp = load(params,'params');
  params = tmp.params;
end

if chr<1 || chr>24, error('chr should be 1-24'); end

% load exome

b=read_table(['/xchip/tcga_scratch/gadgetz/mutation_analysis/' ...
              'whole_exome_refseq_coding.targets.interval_list.no_header'],...
             '%s%f%f%*s%*s',char(9),0);
X.chrn = convert_chr(b.dat{1});
X.start = b.dat{2};
X.end = b.dat{3};

chrlen = load_chrlen;
X2 = reorder_struct(X,X.chrn==chr & X.end<chrlen(chr));

M={};idx=1;
for first=1:1000:slength(X2)
  last=min(slength(X2),first+999);
  fprintf('\n\nTARGETS %d-%d of %d\n\n',first,last,slength(X2));
  try
    result = somcall(t_boomdir,n_boomdir,X2.chrn(first:last),X2.start(first:last),X2.end(first:last),params);
    if ~isempty(result)
      M{idx} = result;
      idx=idx+1;
    end
  catch me
    fprintf('\n\n*****************\n  CHUNK FAILED!\n*****************\n\n');
  end
end

M = cat(1,M{:});
  
if ~exist(outdir,'dir'), mkdir(outdir); end
save([outdir '/chr' num2str(chr) '.mat'],'M');
