function run(cset)

build = 'hg18'

dir = '/xchip/cga1/lawrence/db/mmexpr2';
if ~exist(dir,'dir'), mkdir(dir); end

chrlen = load_chrlen;
all = load_struct([dir '/mm_expr2_genes.txt']);
all = make_numeric(all,'categ');
R = load_refseq(build);
R.chr = convert_chr(R.chr);
R.aidx = match_genelists(upper(R.gene),upper(all.name));

% ADD PROMOTERS TO DATABASE REGIONS
promoter_size = 3000;
R.gene_start = R.tx_start;
R.gene_end = R.tx_end;
for i=1:slength(R)
  switch(R.strand{i})
   case '+', R.gene_start(i) = R.gene_start(i) - promoter_size;
   case '-', R.gene_end(i) = R.gene_end(i) + promoter_size;
  end
end

if ~exist('cset','var'), cset=1:24; end
for csetidx=1:length(cset)
chr = cset(csetidx);
fprintf('\nchr%d\n',chr);
if chr<1 || chr>24, error('chr must be 1-24'); end

if chr==1
  categ_list = [];
  categ_list.num = (0:10)';
  categ_list.name = {'unknown';'poorly_performing_probes';'expression_0';'expression_gt0';...
                     'expression_gt1pct';'expression_gt5pct';'expression_gt25pct';...
                     'expression_gt50pct';'expression_gt75pct';'expression_gt95pct';'expression_gt99pct'};
  save_struct(categ_list,[dir '/categs.txt']);
end

categ = zeros(chrlen(chr),1);
for i=1:slength(R)
  if R.chr(i)~=chr, continue; end
  if isnan(R.aidx(i)), continue; end
  st = R.gene_start(i); en = R.gene_end(i);
  categ(st:en) = max(categ(st:en),all.categ(R.aidx(i)));
end

% SAVE
fname = [dir '/chr' num2str(chr) '.mat'];
save(fname,'categ');
f = fopen(regexprep(fname,'\.mat$','\.txt'),'wt');
fprintf(f,'%d\n',categ);
fclose(f);

end % next chr

