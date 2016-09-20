function run(cset)

if ~exist('cset','var'), cset=1:24; end
for csetidx=1:length(cset)
c = cset(csetidx);
fprintf('\nchr%d\n',c);
if c<1 || c>24, error('c must be 1-24'); end

IGR = 0
INTRON = 1
PROMUTR = 2
EXON = 3

impute_promoters = true;

if imput_promoters
  categ_list.num = [IGR;INTRON;PROMUTR;EXON];
  categ_list.name = {'IGR';'intron';'promoter/UTR';'exon'};
  promoter_size = 3000;
else
  categ_list.num = [IGR;INTRON;PROMUTR;EXON];
  categ_list.name = {'IGR';'intron';'UTR';'exon'};
  promoter_size = 0;
end

if ~exist('loaded_flag','var')
  Rall = load_genedb;
  Rall.chrom = convert_chr(Rall.chrom);
  % ADD PROMOTERS TO DATABASE REGIONS
  Rall.geneStart = Rall.txStart;
  Rall.geneEnd = Rall.txEnd;
  for i=1:slength(Rall)
    switch(Rall.strand{i})
      case '+', Rall.geneStart(i) = Rall.geneStart(i) - promoter_size;
      case '-', Rall.geneEnd(i) = Rall.geneEnd(i) + promoter_size;
    end
  end
  chrlen = load_chrlen;
  loaded_flag = true;
end

ln = chrlen(c);
categ = IGR*ones(ln,1);
R = reorder_struct(Rall,Rall.chrom==c);
for i=1:slength(R)
    if ~mod(i,1000), fprintf('%d/%d ',i,slength(R)); end
    % left UTR/Promoter
    dst = max(1,R.geneStart(i));
    den = min(ln,R.cdsStart(i));
    categ(dst:den) = max(PROMUTR,categ(dst:den));
    % exons and introns
    for e=1:R.exonCount(i)
      est = R.exonStarts{i}(e);
      % intron
      if e>1
        ist = een+1;  % end of last exon
        ien = est-1;  % start of this exon
        dst = max(1,ist); den = min(ln,ien);
        categ(dst:den) = max(INTRON,categ(dst:den));
      end
      % exon
      een = R.exonEnds{i}(e);
      cst = est; cen = een; noncoding = false;
      if R.cdsStart(i)>cen || R.cdsEnd(i)<cst
        noncoding = true;
      elseif R.cdsStart(i)>cst && R.cdsEnd(i)<cen
        cst = R.cdsStart(i); cen = R.cdsEnd(i);
      elseif R.cdsStart(i)>cst; cst = R.cdsStart(i);
      elseif R.cdsEnd(i)<cen; cen = R.cdsEnd(i);
      end
      if ~noncoding
        dst = max(1,cst); den = min(ln,cen);
        categ(dst:den) = max(EXON,categ(dst:den));
      end
    end
    % right UTR/Promoter
    dst = max(1,R.cdsEnd(i));
    den = min(ln,R.geneEnd(i));
    categ(dst:den) = max(PROMUTR,categ(dst:den));
end

if c==1
  fname = '/xchip/tcga_scratch/lawrence/db/zone/categs.txt';
  save_struct(categ_list,fname);
end

fname = ['/xchip/tcga_scratch/lawrence/db/zone/chr' num2str(c) '.mat'];
save(fname,'categ');
f = fopen(regexprep(fname,'\.mat$','\.txt'),'wt');
fprintf(f,'%d\n',categ);
fclose(f);

end

