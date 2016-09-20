function extract_from_wig(samples,targetlist,outname,categfile,build,build_dir,mincateg,maxcateg,wigname,P)
% extract_from_wig(samples,targetlist,outname,categfile,build,build_dir,mincateg,maxcateg,wigname,P)
%
% samples can be a single sample name or a cell array of sample names
%
% targetlist points to a file with one line per target region:  <targname/genename>   <chr>   <start>   <end>
%      can be tab- or space-delimited
%
% outname is the output file that will be created in the sample's home directory
%     ALTERNATE USAGE: outname can be given as a cell array of absolute paths to the output files to be written
% 
%   output file has one line per target:
%     each line will list the number of "covered" bases in the target (based on 8/14 cutoffs) for each category

% categfile (optional) is a category FWB file telling categories to count the coverage in 
%   e.g. /xchip/tcga_scratch/lawrence/db/context/all.fwb
%   if specified, mincateg and maxcateg must also be supplied.
%
% build, e.g. hg19, mm9, canFam2
%
% build_dir, e.g. /xchip/cga/reference/annotation/db/ucsc/hg19/hg19_info.txt
%
% wigname = name of coverage wiggle file to read from (from within the "sample" directory e.g. gbm/0188/wgs)
%     ALTERNATE USAGE: wigname can be given as a cell array of absolute paths to the wiggle files.
%                      some/all members of the cell arrays can be cell arrays of wiggle files to "OR" together.
%             (in this case, collapse_flag cannot be specified as true)
%
% Mike Lawrence 2009-05-12

if nargin>10, error('argument structure has changed: now takes 9 arguments + optional P struct'); end
if exist('P','var') && ~isempty(P) && ~isstruct(P), error('argument structure has changed: 8th argument should be a P struct'); end

if isnumeric(build) || isnumeric(build_dir)
  error('argument structure of extract_from_wig has changed');
end

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'collapse_flag', false);
P = impose_default_value(P,'wait_flag', true);
P = impose_default_value(P,'expect_bin_file',true);
P = impose_default_value(P,'bsub_queue',[]);

demand_file(targetlist);

if ~exist('categfile','var') || isempty(categfile)
  categfile = 'none';
  mincateg = 1;
  maxcateg = 1;
else 
  if ~exist('mincateg','var') || isempty('mincateg') || ~exist('maxcateg','var') || isempty('maxcateg')
    error('If <categfile> is specified, <mincateg> and <maxcateg> must also be specified.\n');
  end
  %% check to see if *.wig is still being used for categfile instead of *.fwb as now required
  if length(categfile)<5, error('<categfile> too short to be a proper filename'); end
  if strcmp(categfile(end-3:end),'.wig')
    oldcategfile = categfile;
    categfile = regexprep(oldcategfile,'\.wig$','.fwb');
    fprintf('Changing %s to %s\n',oldcategfile,categfile);
  end
  demand_file(categfile);
end
if mincateg<1 || maxcateg<1 || mincateg>32767 || maxcateg>32767 || mincateg>maxcateg
  error('Problem with mincateg, maxcateg');
end

if ~exist('wigname','var'), wigname = 'somatic_coverage.wig'; end
%if ~exist('collapse_flag','var'), collapse_flag = false; end
%if ~exist('wait_flag','var'), wait_flag = true; end

if ~iscell(samples), samples = {samples}; end

java_classpath = get_jcp;
java_class = 'org.broadinstitute.cga.tools.seq.ExtractFromWigAllChr';

% DETERMINE INPUT AND OUTPUT FILES

infile = cell(length(samples),1);
outfile = cell(length(samples),1);

% ALTERNATE USAGE MODE (absolute paths, ignoring "samples")
if iscell(outname) || iscell(wigname)
  if ~iscell(outname) || ~iscell(wigname), error('To use absolute-path mode, both outname and wigname must be cell arrays'); end  
  if P.collapse_flag, error('P.collapse_flag is incompatible with absolute-path mode'); end
  if length(outname) ~= length(wigname), error('length(outname)~=length(wigname)'); end
  if length(outname) ~= length(samples), error('length(outname)~=length(samples)'); end
  for s=1:length(samples)
    outfile{s} = outname{s};
    if ischar(wigname{s})
      infile{s} = {wigname{s}};
    elseif iscell(wigname{s})
      infile{s} = wigname{s};
    end
  end
else
  % CLASSICAL USAGE MODE: keyed of "sample" dirs (e.g. gbm/0188/wgs)
  if P.collapse_flag    % collapse patients across centers
    tmp = regexprep(samples,'(.*)/(\d+)/(capture|wgs)-.*','$1/$2/$3');
    [u ui uj] = unique(tmp);
    indir = cell(length(u),1);
    for i=1:length(u), indir{i} = samples(uj==i); end
    samples = tmp(ui);
  else
    indir = cell(length(samples),1);
    for i=1:length(samples), indir{i} = samples(i); end
  end
  bd = '/xchip/cga1/lawrence';
  for s=1:length(samples)
    sample = samples{s};
    outdir = [bd '/' sample];
    if ~exist(outdir,'dir'), fprintf('Creating %s\n',outdir); mkdir(outdir); end
    outfile{s} = [outdir '/' outname];
    infile{s} = cell(length(indir{s}),1);
    for i=1:length(indir{s})
      infile{s}{i} = [bd '/' indir{s}{i} '/' wigname];
    end
  end
end

% SUBMIT JOBS

if P.expect_bin_file
  expected_outfile = regexprep(outfile,'^(.*)$','$1.bin');
else
  expected_outfile = outfile;
end


PP=[];
PP.queue = P.bsub_queue;

did_nothing = true;
all_done = false;
while(~all_done)
  jobs=[];
  cmds = {}; banners = {};
  all_done = true;
  for s=1:length(samples)
    sample = samples{s};
    short_sample = sample_to_short_sample(sample);
    name2 = upper(regexprep(sample,'/','-'));
    name2 = regexprep(name2,'-(ABI|WGS|CAPTURE)$','');
    patient_name = name2;
    % check to see if outfile is older than any infiles
    dout = dir(expected_outfile{s});
    uptodate = true;
    for i=1:length(infile{s})
      demand_file(infile{s}{i});
      din = dir(infile{s}{i});
      if isempty(dout) || dout.datenum<din.datenum, uptodate = false; end
    end
    if uptodate, continue; end
    % create job
    did_nothing = false;
    all_done = false;
    banners{end+1,1} = [short_sample 'XFW'];
    cmds{end+1,1} = ['-R "rusage[mem=7]" "java -Xmx5g -classpath ' java_classpath ' ' java_class ...
        ' ' targetlist ' ' categfile ' ' num2str(mincateg) ' ' num2str(maxcateg)...
        ' ' patient_name ' ' build ' ' build_dir ' ' outfile{s}];
    for i=1:length(infile{s}), cmds{end} = [cmds{end} ' ' infile{s}{i}]; end
    cmds{end}(end+1) = '"';
    % submit jobs as they accumulate
    if length(cmds)>=60
      jobs = [jobs; bsub(cmds,banners,PP)];
      cmds = {}; banners = {};
    end
  end % next sample
  % submit any remaining jobs
  if ~isempty(cmds)
    jobs = [jobs; bsub(cmds,banners,PP)];
    cmds = {}; banners = {};
  end
  if (all_done) break; end
  if (~P.wait_flag) return; end
  bwait(jobs);
end

if did_nothing, fprintf('All files are already up-to-date\n');
else fprintf('All files are now up-to-date\n');
end

