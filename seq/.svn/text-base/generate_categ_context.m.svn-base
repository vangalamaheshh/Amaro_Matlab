function run(cset)

build = 'hg18'

if ~exist('cset','var'), cset=1:24; end
for csetidx=1:length(cset)
c = cset(csetidx);
fprintf('\nchr%d\n',c);

if c<1 || c>24, error('c must be 1-24'); end

categ_list.num = [1;2;3;4];
categ_list.name = {'CpG';'other CG';'AT';'N'};

d = upper(genome_region(c,1,inf),build);
categ = 4*ones(length(d),1);
categ(d=='A' | d=='T') = 3;
categ(d=='C' | d=='G') = 2;
idx = find(d(1:end-1)=='C' & d(2:end)=='G');
categ([idx;idx+1]) = 1;

if c==1
  fname = '/xchip/tcga_scratch/lawrence/db/context/categs.txt';
  save_struct(categ_list,fname);
end

fname = ['/xchip/tcga_scratch/lawrence/db/context/'...
  'chr' num2str(c) '.mat'];
save(fname,'categ');
f = fopen(regexprep(fname,'\.mat$','\.txt'),'wt');
fprintf(f,'%d\n',categ);
fclose(f);

end
