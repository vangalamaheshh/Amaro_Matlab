function run(cset)

dir = '/xchip/cga1/lawrence/db/islands';
if ~exist(dir,'dir'), mkdir(dir); end

chrlen = load_chrlen;
I = load_struct([dir '/CpG_islands.txt']); % downloaded from UCSC ("CpG islands track")
I.chr = convert_chr(I.chrom);
I = reorder_struct(I,~isnan(I.chr));
I.start = str2double(I.chromStart);
I.end = str2double(I.chromEnd);

shoresize = 2000;

if ~exist('cset','var'), cset=1:24; end
for csetidx=1:length(cset)
  chr = cset(csetidx);
  fprintf('\nchr%d\n',chr);
  if chr<1 || chr>24, error('chr must be 1-24'); end

  if chr==1
    categ_list = [];
    categ_list.num = (1:3)';
    categ_list.name = {'sea';'shore';'island'};
    save_struct(categ_list,[dir '/categs.txt']);
  end

  len = chrlen(chr);  
  categ = ones(len,1);
  Ic = reorder_struct(I,I.chr==chr);
  for i=1:slength(Ic)
    st1 = max(1,min(len,Ic.start(i)));
    en1 = max(1,min(len,Ic.end(i)));
    st2 = max(1,min(len,Ic.start(i)-shoresize));
    en2 = max(1,min(len,Ic.end(i)+shoresize));
    categ(st2:en2) = max(2,categ(st2:en2));
    categ(st1:en1) = max(3,categ(st1:en1));
  end
        
  % SAVE
  fname = [dir '/chr' num2str(chr) '.mat'];
  save(fname,'categ');
  f = fopen(regexprep(fname,'\.mat$','\.txt'),'wt');
  fprintf(f,'%d\n',categ);
  fclose(f);

end % next chr
