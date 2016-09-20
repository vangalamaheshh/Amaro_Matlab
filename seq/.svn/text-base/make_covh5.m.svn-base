function make_covh5(samples)
% make_covh5(samples)
%
% runs igv tool "coverage.sh" to create coverage *.h5 file for the tumor and normal bams specified
%
% Mike Lawrence 2009-08-18

if ~exist('samples','var'), error('<samples> is required'); end
if ~iscell(samples), samples = {samples}; end

if ~exist('P','var'), P=[]; end
if ~isempty(P) && ~isstruct(P), error('Second parameter should be a P struct'); end

P = impose_default_value(P,'cancer_sample',true);
P = impose_default_value(P,'coverage_window_size',25);

scriptname = '/xchip/tcga_scratch/ng/tools/igv/coverage.sh';

try

if P.cancer_sample, tn = {'normal','tumor'}; else tn = {'sample'}; end

did_nothing = true;
all_done = false;
while(~all_done)
  jobs=[]; all_done = true;
  for i=1:length(samples)
    sample = samples{i};
    direc = ['/xchip/tcga_scratch/lawrence/' sample];
    for t=1:length(tn)
      stem = [direc '/' tn{t}];
      bam = [stem '.bam'];
      cov = [bam '.h5'];
      dbam = dir(bam);
      if isempty(dbam)
        fprintf('ERROR: %s not found!\n',bam);
        continue
      end
      dcov = dir(cov);
      if ~isempty(dcov)
        if datenum(dcov.date)>=datenum(dbam.date)
          continue;          % cov is already up-to-date
      end,end
      all_done = false;
      did_nothing = false;
      cmd = ['"' scriptname ' ' bam ' ' direc ' hg18 -w ' num2str(P.coverage_window_size) '"'];
      banner = [sample_to_short_sample(sample) 'COVH5' tn{t}(1)];
      jobs = [jobs; bsub(cmd,banner)];
  end,end
  if (~all_done), bwait(jobs); else break; end
end

if did_nothing, fprintf('Coverage H5''s are already up-to-date for all those samples.\n'); end

catch me, excuse(me); end

