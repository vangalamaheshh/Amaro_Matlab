function extract_from_cbb(samples,targetlist,outname,categdir,numcategs)
% extract_from_cbb(samples,targetlist,outname,[categdir,numcategs])
%
% samples can be a single sample name or a cell array of sample names
%
% targetlist points to a file with one line per target region:  <targname/genename>   <chr>   <start>   <end>
%      can be tab- or space-delimited
%
% outname is the output file that will be created in the sample's home directory
% 
% output file has one line per target:
%   each line lists the number of "covered" bases in the target (based on 8/14 cutoffs)
%      for each category
%
% categdir (optional) is a directory containing category files
%   (one file per chromosome, one line per basepair, telling which category to count the coverage in) 
%   e.g. /xchip/tcga_scratch/lawrence/db/context
%
%   if specified, numcategs must also be supplied.
%
% Mike Lawrence 2009-05-12

if ~exist('categdir','var') | isempty('categdir')
  categdir = 'none';
  numcategs = 1; 
else 
  if ~exist('numcategs','var') | isempty('numcategs')
    error('If <categdir> is specified, <numcategs> must also be specified.\n');
  end
end

if ~iscell(samples), samples = {samples}; end

script = '/home/radon00/lawrence/CancerGenomeAnalysis/trunk/matlab/seq/extract_from_cbb.pl';

try

% first see if final output file exists for any of these samples
already_done = {};
for s = 1:length(samples)
  sample = samples{s};
  direc = ['/xchip/cga1/lawrence/' sample];
  tumdir = [direc '/tumor_cbb'];
  normdir = [direc '/normal_cbb'];
  outfile = [direc '/' outname];
  dout = dir(outfile);
  if ~isempty(dout)
    uptodate = true;
    for c=1:24
      tumfile = [tumdir '/chr' num2str(c) '.cbb'];
      normfile = [normdir '/chr' num2str(c) '.cbb'];
      dtum = dir(tumfile);
      if isempty(dtum) | dtum.bytes<50000000, fprintf('ERROR: %s incomplete or missing\n', tumfile); uptodate=false;
      else if dtum.datenum>dout.datenum, uptodate = false; end; end
      dnorm = dir(normfile);
      if isempty(dnorm) | dnorm.bytes<50000000, fprintf('ERROR: %s incomplete or missing\n', normfile); uptodate=false;
      else if dnorm.datenum>dout.datenum, uptodate = false; end; end
    end
    if uptodate, already_done = [already_done; sample]; end
  end
end
samples = setdiff(samples,already_done);
if length(samples)==0, fprintf('All files are up-to-date\n'); return; end

% submit per-chromosome jobs for each sample that needs processing

all_done = false;
while(~all_done)
  jobs=[];
  all_done = true;
  for s = 1:length(samples)
    sample = samples{s};
    short_sample = sample_to_short_sample(sample);
    direc = ['/xchip/cga1/lawrence/' sample];
    tumdir = [direc '/tumor_cbb'];
    normdir = [direc '/normal_cbb'];
    outfile = [direc '/' outname];
    outstem = [outfile '_tmp_'];
    for c=1:24
      outfile = [outstem '_chr' num2str(c)];
      tumfile = [tumdir '/chr' num2str(c) '.cbb'];
      normfile = [normdir '/chr' num2str(c) '.cbb'];
      dtum = dir(tumfile);
      if isempty(dtum) | dtum.bytes<50000000, fprintf('ERROR: %s incomplete or missing\n', tumfile); end
      dnorm = dir(normfile);
      if isempty(dnorm) | dnorm.bytes<50000000, fprintf('ERROR: %s incomplete or missing\n', normfile); end
      dout = dir(outfile);
      if ~isempty(dout)
        if dout.datenum>=dnorm.datenum & dout.datenum>=dtum.datenum
          continue;   % file already exists and is up-to-date
        else
          fprintf('%s needs to be refreshed:\n',outfile);
        end
      end
      all_done = false;
      fprintf('Creating %s\n',outfile);
      banner = [short_sample 'XFC' num2str(c)];
      cmd = ['perl ' script ...
        ' ' targetlist ' ' tumdir ' ' normdir ' ' categdir ' ' num2str(numcategs) ' ' outstem ' ' num2str(c)];
%      jobs = [jobs;bsub(['-R "rusage[mem=8]" "' cmd '"'],banner)];
 jobs = [jobs;bsub(['"' cmd '"'],banner)];

    end % next chr
  end % next sample
  if (all_done) break; end
  bwait(jobs);
end

% concatenate and delete temporary files

fprintf('Concatenate and delete temporary files:\n');
for s = 1:length(samples)
    sample = samples{s};
    fprintf('  Samples %s\n',sample);
    direc = ['/xchip/cga1/lawrence/' sample];
    outfile = [direc '/' outname];
    tmpfile = [direc '/X' outname];
    outstem = [outfile '_tmp_'];

    system(['cat ' outstem '_chr* > ' tmpfile]);
    if exist(tmpfile,'file')
      system(['rm ' outstem '_chr*']);
      system(['mv ' tmpfile ' ' outfile]);
    else
      fprintf('Problem concatenating files for sample %s\n',sample);
      keyboard
    end
end

catch me, excuse(me); end
