function extract_from_fwb(samples,targetlist,outname,categfile,mincateg,maxcateg,fwbname,collapse_flag,wait_flag)
% extract_from_fwb(samples,targetlist,outname,[categfile,mincateg,maxcateg,fwbname,collapse_flag,wait_flag])
%
% samples can be a single sample name or a cell array of sample names
%
% targetlist points to a file with one line per target region:  <targname/genename>   <chr>   <start>   <end>
%      can be tab- or space-delimited
%
%   outname must be given as a cell array of absolute paths to the output files to be written
% 
%   output file has one line per target:
%     each line will list the number of "covered" bases in the target (based on 8/14 cutoffs) for each category

%   categfile (optional) is a category fwb file telling categories to count the coverage in 
%     e.g. /xchip/tcga_scratch/lawrence/db/context/all.fwb
%     if specified, mincateg, maxcateg, and categwidth must also be supplied.
%
%   fwbname must be given as a cell array of absolute paths to the fwb files.
%                      some/all members of the cell arrays can be cell arrays of fwb files to "OR" together.
%                      (in this case, collapse_flag cannot be specified as true)
%
% Mike Lawrence 2010-10-13

demand_file(targetlist);

if ~exist('categfile','var') || isempty('categfile')
  categfile = 'none';
  mincateg = 1;
  maxcateg = 1;
else 
  if ~exist('mincateg','var') || isempty('mincateg') || ~exist('maxcateg','var') || isempty('maxcateg')
    error('If <categfile> is specified, <mincateg> and <maxcateg> must also be specified.\n');
  end
  if isempty(grepi('.fwb$',categfile,1)), error('categfile needs to be an FWB'); end
  demand_file(categfile);
end
if mincateg<1 || maxcateg<1 || mincateg>maxcateg
  error('Problem with mincateg, maxcateg');
end

if ~exist('fwbname','var'), error('need fwbname'); end
if ~exist('collapse_flag','var'), collapse_flag = false; end
if ~exist('wait_flag','var'), wait_flag = true; end

if ~iscell(samples), samples = {samples}; end

java_classpath = get_jcp;
java_class = 'org.broadinstitute.cga.tools.seq.ExtractFromFwbAllChr';

% DETERMINE INPUT AND OUTPUT FILES

infile = cell(length(samples),1);
outfile = cell(length(samples),1);

if iscell(outname) || iscell(fwbname)
  if ~iscell(outname) || ~iscell(fwbname), error('To use absolute-path mode, both outname and fwbname must be cell arrays'); end  
  if collapse_flag, error('collapse_flag is incompatible with absolute-path mode'); end
  if length(outname) ~= length(fwbname), error('length(outname)~=length(fwbname)'); end
  if length(outname) ~= length(samples), error('length(outname)~=length(samples)'); end
  for s=1:length(samples)
    outfile{s} = outname{s};
    if ischar(fwbname{s})
      infile{s} = {fwbname{s}};
    elseif iscell(fwbname{s})
      infile{s} = fwbname{s};
    end
  end
else
  error('classical usage mode no longer supported');
end

% SUBMIT JOBS

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
    dout = dir(outfile{s});
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
    banners{end+1,1} = [short_sample 'XFF'];
    cmds{end+1,1} = ['"java -classpath ' java_classpath ' ' java_class ...
        ' ' targetlist ' ' categfile ' ' num2str(mincateg) ' ' num2str(maxcateg) ...
        ' ' patient_name ' ' outfile{s}];
    for i=1:length(infile{s})
      if isempty(grepi('.fwb',infile{s}{i},1)), error('covfiles need to be FWB'); end
      cmds{end} = [cmds{end} ' ' infile{s}{i}];
    end
    cmds{end}(end+1) = '"';
    % submit jobs as they accumulate
    if length(cmds)>=60
      jobs = [jobs; bsub(cmds,banners)];
      cmds = {}; banners = {};
    end
  end % next sample
  % submit any remaining jobs
  if ~isempty(cmds)
    jobs = [jobs; bsub(cmds,banners)];
    cmds = {}; banners = {};
  end
  if (all_done) break; end
  if (~wait_flag) return; end
  bwait(jobs);
end

if did_nothing, fprintf('All files are already up-to-date\n');
else fprintf('All files are now up-to-date\n');
end

