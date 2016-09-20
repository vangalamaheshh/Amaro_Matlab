function get_mutrates(sample)
% get_mutrates(sample)
%
% sample
%    e.g.    gbm/0188
%            cll/je
%            ov/0725
%            ov/0751

try

readlength = getBAMFileReadLength(sample);

for i=1:2
  if i==1, tn='normal';else tn='tumor'; end
  cbbDir{i} = ['/xchip/tcga_scratch/lawrence/' sample '/' tn '_cbb'];
end

% load mutation data

name2 = regexprep(upper(sample),'/','-');
fname = ['/xchip/tcga_scratch/ng/' name2 '/wgs/mut/mutation_reports.maf.annotated'];
M = load_struct(fname);
M.cc = convert_chr(M.chr);
M.start = str2double(M.start);
M.end = str2double(M.end);


% make sure all needed files exist

use=true(24,1);
for chr=1:24
  file = ['chr' num2str(chr) '.cbb'];
  for i=1:2
    f = [cbbDir{i} '/' file];
    if ~exist(f), fprintf('%s does not exist\n',f); use(chr)=false;end
    d = dir(f);
    if d.bytes<50000000, fprintf('%s is too short\n',f);use(chr)=false; end
  end
end

% load reference data

R = load_genedb;
R.chr = convert_chr(R.chrom);
% note: USCS zero-based/one-based/half-closed set convention
%   is already corrected throughout to inclusive one-based convention

for chr=1:24, if ~use(chr), continue; end
  fprintf('Analyzing chromosome %d\n',chr);

  % load genome
  fprintf('Loading reference genomic sequence\n');
  G = genome_region(chr,1,inf);
  clen = length(G);

  % partition genome into contexts 11=A_A 12=A_C 13=A_G 14=A_T 21=A_C ... 44=T_T
  C = upper(G);
  C(C=='A')=1; C(C=='C')=2; C(C=='G')=3; C(C=='T')=4; C(C=='N')=0;
  C = double(C); C(C==0) = NaN;

  % partition genome into transcribed/untrabscribed and intron/exon regions
  % 0=IGR  1=intron  2=UTR  3=exon
  D = zeros(clen,1);
  idx = find(R.chr==chr);
  fprintf('Parsing gene ');
  for j=1:length(idx), i=idx(j);
    if ~mod(j,100) fprintf('%d/%d ',j,length(idx)); end
    % left UTR
    dst = min(1,R.txStart(i));
    den = min(clen,R.cdsStart(i));
    D(dst:den) = max(2,D(dst:den));
    % exons and introns
    for e=1:R.exonCount(i)
      est = R.exonStarts{i}(e);
      % intron
      if e>1
        ist = max(1,een+1);  % end of last exon
        ien = min(clen,est-1);  % start of this exon
        D(ist:ien) = max(1,D(ist:ien));
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
        dst = max(1,cst); den = min(clen,cen);
        D(dst:den) = max(3,D(dst:den));
      end
    end
    % right UTR
    dst = max(1,R.cdsEnd(i));
    den = min(clen,R.txEnd(i));
    D(dst:den) = max(2,D(dst:den));
  end
  fprintf('Done.\n');

  % retrieve mutations
  Mc = reorder_struct(M,M.cc==chr);
  Mc.region_type = D(Mc.start);

  % load CoverageByBase data
  fprintf('Loading CoverageByBase data\n');
  file = ['chr' num2str(chr) '.cbb'];
  N = sscanf(load_textfile([cbbDir{1} '/' file]),'%d');
  T = sscanf(load_textfile([cbbDir{2} '/' file]),'%d');

  T(clen)=0;
  N(clen)=0;
  T=T(1:clen);
  N=N(1:clen);

  fprintf('Sample %s\n',sample);

  for kind=0:4   % 0=IGR  1=intron  2=UTR  3=exon   4=total

    if kind<4
      idxx = find(D==kind);
      nmuts(chr,kind+1) = sum(Mc.region_type==kind);
    else
      idxx=1:length(D);
      nmuts(chr,kind+1) = slength(Mc);
    end

    tumcov(chr,kind+1) = sum(N(idxx));
    normcov(chr,kind+1) = sum(T(D==idxx));

    if length(T)~=length(N)
      fprintf('length disagree\n');keyboard
    else
      callcov(chr,kind+1) = sum(T(idxx)>=14 & N(idxx)>=8);
    end
  end
end

for kind=0:4
  switch kind
    case 0: fprintf('IGRs\n');
    case 1: fprintf('introns\n');
    case 2: fprintf('UTRs/promoters\n');
    case 3: fprintf('exons\n');
    case 4: fprintf('Total\n');
  end

  tsum = sum(tumcov(use,kind+1));
  nsum = sum(normcov(use,kind+1));
  csum = sum(callcov(use,kind+1));
  totlen = sum(chrlen(use));

  treads = tsum / readlength * (sum(chrlen)/totlen);
  nreads = nsum / readlength * (sum(chrlen)/totlen);

  tdepth = tsum / totlen;
  ndepth = nsum / totlen;
  avgc = csum / totlen;

  nm = sum(nmuts,kind);
  mutrate = nm / csum;

  fprintf('  Total reads: %d tumor,  %d normal\n',treads,nreads);
  fprintf('  Average read depth:  %0.2f tumor,  %0.2f normal\n',tdepth,ndepth);
  fprintf('  Fraction of genome callable for somatic mutations:  %0.2f\n',avgc);
  fprintf('  Mutation rate %0.2f\n',mutrate);

catch me, excuse(me); end

