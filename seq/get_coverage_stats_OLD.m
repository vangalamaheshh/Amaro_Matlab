function get_coverage_stats(sample,filters)
% get_coverage_stats(sample,filters)
%
% sample
%    e.g.    gbm/0188/wgs
%            cll/je/wgs
%            ov/0725/c2k
%            ov/0751/c6k
%    (multiple samples not currently supported)
%
% filters
%    e.g.    context
%            c2k 
%            c6k
%    (can be a cell array of strings)

try

max_ncat = 100;
if ~exist('filters','var'), filters=[]; end
if ~isempty(filters)
  if ~iscell(filters), filters = {filters}; end
  nf = length(filters);
  F = cell(nf,1);
  for f=1:nf
    F{f}.categ_list = load_struct(['/xchip/tcga_scratch/lawrence/db/' filters{f} '/categs.txt'],'%f%s');
    F{f}.ncat = slength(F{f}.categ_list);
    if F{f}.ncat>max_ncat, error('Filter %s has more than max=%d categories.',filters{f},max_ncat); end
  end
else
  nf = 0;
end

readlength = getBAMFileReadLength(sample);

for i=1:2
  if i==1, tn='normal';else tn='tumor'; end
  cbbDir{i} = ['/xchip/tcga_scratch/lawrence/' sample '/' tn '_cbb'];
end

% make sure all needed files exist

use=true(24,1);
for chr=1:24
  file = ['chr' num2str(chr) '.cbb'];
  for i=1:2
    f = [cbbDir{i} '/' file];
    if ~exist(f), fprintf('WARNING: %s does not exist\n',f); use(chr)=false;end
    d = dir(f);
    if d.bytes<50000000, fprintf('WARNING: %s is too short\n',f);use(chr)=false; end
  end
end

if ~any(use), error('No coverage files have been generated!'); end
if ~all(use), error('Not all coverage files have been generated!'); end

chrlen = load_chrlen;

fprintf('Analyzing chromosomes: \n');

nz = nf+1;
tumcov = nan(24,nz,max_ncat);
normcov = nan(24,nz,max_ncat);
callcov = nan(24,nz,max_ncat);
territory = nan(24,nz,max_ncat);

for chr=1:24, if ~use(chr), continue; end
  fprintf('%d ',chr);
  file = ['chr' num2str(chr) '.cbb'];

  N = sscanf(load_textfile([cbbDir{1} '/' file]),'%d');
  T = sscanf(load_textfile([cbbDir{2} '/' file]),'%d');

  m = chrlen(chr);
  T(m+1)=0; N(m+1)=0;
  T=T(1:m); N=N(1:m);

  % unfiltered
  tumcov(chr,1,1) = sum(T);
  normcov(chr,1,1) = sum(N);
  territory(chr,1,1) = m;
  callcov(chr,1,1) = sum(T>=14 & N>=8);

  % filters
  if nf>0
    for f=1:nf
      V = load(['/xchip/tcga_scratch/lawrence/db/' filters{f} '/chr' num2str(chr)]);
      for k=1:F{f}.ncat
        mask = (V.categ==F{f}.categ_list.num(k));
        mask(m+1)=0; mask=mask(1:m);
        tumcov(chr,f+1,k) = sum(T(mask));
        normcov(chr,f+1,k) = sum(N(mask));
        territory(chr,f+1,k) = sum(mask);
        callcov(chr,f+1,k) = sum(T(mask)>=14 & N(mask)>=8);
  end,end,end
end

% REPORT

out = fopen(['/xchip/tcga_scratch/lawrence/' sample '/coverage.stats'],'wt');
fprintf(out,'Sample %s\n',sample);
for z=1:nz

  if z==1, fprintf(out,'  WHOLE GENOME\n'); ncat = 1;
  else fprintf(out,'  FILTERED BY %s\n',upper(filters{z-1})); ncat = F{z-1}.ncat; end

  for k=1:ncat
    if z>1
      fprintf(out,'    Category %d = %s\n',F{z-1}.categ_list.num(k),F{z-1}.categ_list.name{k});
    end
    tsum = fullsum(tumcov(use,z,k));
    nsum = fullsum(normcov(use,z,k));
    csum = fullsum(callcov(use,z,k));
    terr = fullsum(territory(use,z,k));

    treads = round(tsum / readlength);
    nreads = round(nsum / readlength);

    tdepth = tsum / terr;
    ndepth = nsum / terr;
    avgc = csum / terr;

    fprintf(out,'      Territory (basepairs):  %.0f\n',terr);
    fprintf(out,'      Number of sequenced bases:  %.0f tumor,  %.0f normal\n',tsum,nsum);
    fprintf(out,'      Number of reads:  %.0f tumor,  %.0f normal\n',treads,nreads);
    fprintf(out,'      Average read depth:  %0.2f tumor,  %0.2f normal\n',tdepth,ndepth);
    fprintf(out,'      Fraction of territory callable for somatic mutations:  %0.2f\n',avgc);
    fprintf(out,'\n');
  end
end
fclose(out);

catch me, excuse(me); end
