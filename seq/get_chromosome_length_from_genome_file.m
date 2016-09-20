function filesize = get_chromosome_length_from_genome_file(cset,build)
%
% build can be hg17, hg18, hg19, mm9, OR absolute path to directory of chr*.txt files

if ~exist('build', 'var') || ~ischar(build)
  error('Must provide genome build.  Can be hg17, hg18, hg19, mm9, OR absolute path to directory of chr*.txt files')
end

if build(1) ~= '/'
  dirname = ['/xchip/cga1/annotation/db/ucsc/' build '/'];
  %dirname = ['/xchip/tcga/gbm/analysis/lawrence/genome/' build '/'];
else
  dirname = build;
  if dirname(end)~='/', dirname = [dirname '/']; end
  clear build;
end

filesize = nan(length(cset),1);
for cseti=1:length(cset)
c = cset(cseti);

if isnumeric(c)
  if length(c)>1, error('multiple chromosomes not supported'); end
  if c==23, c = 'X';
  elseif c==24, c = 'Y';
  elseif c==25, c = 'M';
  else c = num2str(c); end
end
if iscell(c) && length(c)==1, c=c{1}; end
if ~ischar(c), error('chromsome should be numeric or a string'); end
if ~strncmp(c,'chr',3), c = ['chr' c]; end
if strcmp(c,'chr23'), c = 'chrX'; end
if strcmp(c,'chr24'), c = 'chrY'; end

ok = {'chr0';'chr10_random';'chr10';'chr11_random';'chr11';'chr12';'chr13_random'; ...
     'chr13';'chr14';'chr15_random';'chr15';'chr16_random';'chr16';'chr17_random'; ...
      'chr17';'chr18_random';'chr18';'chr19_random';'chr19';'chr1_random';'chr1'; ...
      'chr20';'chr21_random';'chr21';'chr22_h2_hap1';'chr22_random';'chr22';'chr23'; ...
      'chr24';'chr2_random';'chr2';'chr3_random';'chr3';'chr4_random';'chr4'; ...
      'chr5_h2_hap1';'chr5_random';'chr5';'chr6_cox_hap1';'chr6_qbl_hap2';'chr6_random'; ...
      'chr6';'chr7_random';'chr7';'chr8_random';'chr8';'chr9_random';'chr9';'chrM'; ...
      'chrX_random';'chrX';'chrY'};

if ~ismember(c,ok), error('%s.txt is not a valid genome file',c); end

fname = [dirname c '.txt'];

for attempt=1:100,try

  d = dir(fname);
  if length(d)~=1, error('Couldn''t find genome file %s',fname); end
  filesize (cseti)= d.bytes;
  break

  catch me
    disp(me); disp(me.message);
    fprintf('Waiting ten seconds and trying again...\n');
    pause(10);   % wait ten seconds and try again
  end

end   % next try

end % next cseti (next chr)
    
    
