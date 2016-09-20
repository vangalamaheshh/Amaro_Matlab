function get_targcovpct(samples,targdir,outfilestem,threshold,chrset)
% get_targcovpct(samples,targdir,outfilestem,[chr])
%
% targdir = directory of files defining the targetset
%    e.g.    c2k 
%            c6k
%
% outfilestem
%    e.g.    targcovpct.txt
%
%    --> will write:    tumor_targcovpct.txt
%                      normal_targcovpct.txt
%
% threshold = minimum # of reads to count as 'covered'
%    e.g.    20
%
% chrset = chromosome set to count
%    e.g.    1:23
%
%    default = 1:24
%
% Mike Lawrence 2009-08-13

if ~exist('targdir','var'), error('need targdir'); end
if ~exist('outfilestem','var'), error('need outfilestem'); end
if ~exist('threshold','var'), error('need threshold'); end
if ~exist('chrset','var'), chrset = 1:24; end

try

if ~iscell(samples), samples = {samples}; end

scriptname = '/home/radon00/lawrence/CancerGenomeAnalysis/trunk/matlab/seq/targcovpct.pl';
basedir = '/xchip/tcga_scratch/lawrence';

tn = {'tumor','normal'};

for i=1:length(tn)
  cbbdir{i} = [tn{i} '_cbb'];
  outfile{i} = [tn{i} '_' outfilestem];
  tempdir{i} = [outfile{i} '_temp'];
end
targdir = [basedir '/db/' targdir];

% see if any final files are finished

already_done = false(length(samples),length(tn));
for s=1:length(samples)
  sample = samples{s};
  direc = [basedir '/' sample];
  for i=1:length(tn)
    dout = dir([direc '/' outfile{i}]);
    if ~isempty(dout) && dout.bytes>0
      already_done(s,i) = true;
      for cidx=1:length(chrset)
        c = chrset(cidx);
        cbbfile = [direc '/' cbbdir{i} '/chr' num2str(c) '.cbb'];
        dcbb = dir(cbbfile);
        if isempty(dcbb) || dcbb.datenum>dout.datenum, already_done(s,i) = false; break; end
end,end,end,end

if all(already_done(:)), fprintf('All files are up-to-date\n'); return; end

% make sure all need files exist

skip = false(length(samples),1);

for s=1:length(samples)
  sample = samples{s};
  direc = [basedir '/' sample];
  for i=1:length(tn)
    if already_done(s,i), continue; end
    for cidx=1:length(chrset)
      c = chrset(cidx);
      cbbfile = [direc '/' cbbdir{i} '/chr' num2str(c) '.cbb'];
      dcbb = dir(cbbfile);
      if isempty(dcbb) || isempty(dcbb.bytes) || dcbb.bytes==0
        fprintf('%s: file not found or empty\n',cbbfile);
        skip(samples) = true; break;   % skip this sample
end,end,end,end

for cidx=1:length(chrset)
  c = chrset(cidx);
  targfile = [targdir '/chr' num2str(c) '.txt'];
  d = dir(targfile);
  if isempty(d) || isempty(d.bytes) || d.bytes==0
    error('%s: file not found or empty\n',targfile);  % can't continue
end,end

% send out all jobs

all_done = false;
while(~all_done)
  jobs = [];
  all_done = true;
  for s=1:length(samples)
    if skip(s), continue; end
    sample = samples{s};
    direc = [basedir '/' sample];
    ss = sample_to_short_sample(sample);
    for i=1:length(tn)
      if already_done(s,i), continue; end
      thistempdir = [direc '/' tempdir{i}];
      if ~exist(thistempdir,'dir'), mkdir(thistempdir); end
      for cidx=1:length(chrset)
        c = chrset(cidx);
        thiscbbfile = [direc '/' cbbdir{i} '/chr' num2str(c) '.cbb'];
        thisoutfile = [thistempdir '/chr' num2str(c) '.txt'];
        thistargfile = [targdir '/chr' num2str(c) '.txt'];
        if exist(thisoutfile,'file')
          dout = dir(thisoutfile);
          if ~isempty(dout) && ~isempty(dout.bytes) && dout.bytes>0
            dcbb = dir(thiscbbfile);
            if dout.datenum>=dcbb.datenum, continue;
            else fprintf('%s needs to be refreshed:\n',thisoutfile); end
          end
        end
        all_done = false;
        banner = [ss 'TCP' tn{i}(1) num2str(c)];
        cmd = ['"perl ' scriptname ' ' thistargfile ' ' thiscbbfile ' ' num2str(threshold) ' ' thisoutfile '"'];
        jobs = [jobs; bsub(cmd,banner)];
      end % next chr
    end % next tn
  end % next sample
  if (~all_done) bwait(jobs); else break; end
end

% integrate statistics across chromosomes and write final files
for s=1:length(samples)
  if skip(s), continue; end
  sample = samples{s};
  direc = [basedir '/' sample];
  for i=1:length(tn)
    if already_done(s,i), continue; end
    thistempdir = [direc '/' tempdir{i}];
    S = cell(length(chrset),1);
    targ = 0;
    targcov = 0;
    for cidx=1:length(chrset)
      c = chrset(cidx);
      thisoutfile = [thistempdir '/chr' num2str(c) '.txt'];
      x = load_textfile(thisoutfile);
      tmp = sscanf(x,'%f');
      if length(tmp)<2, error('weirdness processing %s!',thisoutfile); end
      targ = targ + tmp(1);
      targcov = targcov + tmp(2);
    end
    save_textfile(sprintf('%d\t%d\t%.2f%%\n',targ,targcov,100*targcov/targ),...
      [direc '/' outfile{i}]);
  end % next tn
end % next sample

fprintf('Finished updating targcovpct files\n');

catch me, excuse(me); end
