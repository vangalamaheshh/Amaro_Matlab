function get_coverage_by_gene(sample)
% get_coverage_by_gene(sample)
%
% for all genes in Refseq, returns number of covered bases in exon regions +/- 2 bp
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

if ~any(use), error('No coverage files have been generated!'); end

chrlen = load_chrlen;

fprintf('Analyzing chromosomes: \n');
for chr=1:24, if ~use(chr), continue; end
  fprintf('%d ',chr);
  file = ['chr' num2str(chr) '.cbb'];

  N = sscanf(load_textfile([cbbDir{1} '/' file]),'%d');
  T = sscanf(load_textfile([cbbDir{2} '/' file]),'%d');

  m = chrlen(chr);
  T(m)=0;
  N(m)=0;
  T=T(1:m);
  N=N(1:m);
  tumcov(chr,1) = sum(N);
  normcov(chr,1) = sum(T);

  if length(T)~=length(N)
    fprintf('length disagree\n');keyboard
  else
    callcov(chr,1) = sum(T>=14 & N>=8);
  end
end

tsum = sum(tumcov(use));
nsum = sum(normcov(use));
csum = sum(callcov(use));
totlen = sum(chrlen(use));

treads = tsum / readlength * (sum(chrlen)/totlen);
nreads = nsum / readlength * (sum(chrlen)/totlen);

tdepth = tsum / totlen;
ndepth = nsum / totlen;
avgc = csum / totlen;

fprintf('Sample %s\n',sample);
fprintf('Total reads: %d tumor,  %d normal\n',treads,nreads);
fprintf('Average read depth:  %0.2f tumor,  %0.2f normal\n',tdepth,ndepth);
fprintf('Fraction of genome callable for somatic mutations:  %0.2f\n',avgc);

catch me, excuse(me); end
