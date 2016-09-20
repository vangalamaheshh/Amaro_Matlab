function run(cset)

if ~exist('cset','var'), cset=1:24; end
for csetidx=1:length(cset)
c = cset(csetidx);
fprintf('\nchr%d\n',c);
if c<1 || c>24, error('c must be 1-24'); end

IGR = 0;
INTRON = 1;
PROMUTR = 2;
EXON = 3;

NEITHER = 0;
PLUS = 1;
MINUS = 2;
BOTH = 3;

impute_promoters = true;

if impute_promoters
  zone_list.num = [IGR;INTRON;PROMUTR;EXON];
  zone_list.name = {'IGR';'intron';'promoter/UTR';'exon'};
  promoter_size = 3000;
else
  zone_list.num = [IGR;INTRON;PROMUTR;EXON];
  zone_list.name = {'IGR';'intron';'UTR';'exon'};
  promoter_size = 0;
end

strand_list.num = [NEITHER;PLUS;MINUS;BOTH];
strand_list.name = {'neither_transcribed';'(+)transcribed';'(-)transcribed';'both_transcribed'};

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
zone = IGR*ones(ln,1);
strand = NEITHER*ones(ln,1);
R = reorder_struct(Rall,Rall.chrom==c);
for i=1:slength(R)
    if ~mod(i,1000), fprintf('%d/%d ',i,slength(R)); end
    % strand transcribed
    dst = max(1,R.txStart(i));
    den = min(ln,R.txEnd(i));
    switch(R.strand{i})
      case '+'
        strand(dst:den) = bitor(PLUS,strand(dst:den));
      case '-'
        strand(dst:den) = bitor(MINUS,strand(dst:den));
    end
    % left UTR/Promoter
    if R.cdsStart(i)<R.cdsEnd(i)
      dst = max(1,R.geneStart(i));
      den = min(ln,R.cdsStart(i));
    else   % noncoding
      dst = max(1,R.geneStart(i));
      den = min(ln,R.cdsEnd(i));
    end
    zone(dst:den) = max(PROMUTR,zone(dst:den));
    % exons and introns
    for e=1:R.exonCount(i)
      est = R.exonStarts{i}(e);
      % intron
      if e>1
        ist = een+1;  % end of last exon
        ien = est-1;  % start of this exon
        dst = max(1,ist); den = min(ln,ien);
        zone(dst:den) = max(INTRON,zone(dst:den));
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
        zone(dst:den) = max(EXON,zone(dst:den));
      end
    end
    % right UTR/Promoter
    dst = max(1,R.cdsEnd(i));
    den = min(ln,R.geneEnd(i));
    zone(dst:den) = max(PROMUTR,zone(dst:den));
end

if ~exist('/xchip/tcga_scratch/lawrence/db/zone','dir'), mkdir('/xchip/tcga_scratch/lawrence/db/zone'); end
if ~exist('/xchip/tcga_scratch/lawrence/db/strand','dir'), mkdir('/xchip/tcga_scratch/lawrence/db/strand'); end

if c==1
  save_struct(zone_list,'/xchip/tcga_scratch/lawrence/db/zone/categs.txt');
  save_struct(strand_list,'/xchip/tcga_scratch/lawrence/db/strand/categs.txt');
end

fname = ['/xchip/tcga_scratch/lawrence/db/zone/chr' num2str(c) '.mat'];
save(fname,'zone');
f = fopen(regexprep(fname,'\.mat$','\.txt'),'wt');
fprintf(f,'%d\n',zone);
fclose(f);

fname = ['/xchip/tcga_scratch/lawrence/db/strand/chr' num2str(c) '.mat'];
save(fname,'strand');
f = fopen(regexprep(fname,'\.mat$','\.txt'),'wt');
fprintf(f,'%d\n',strand);
fclose(f);

end

