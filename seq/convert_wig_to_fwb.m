function convert_wig_to_fwb(targetlist,outname,wigname)
% convert_wig_to_fwb(targetlist,outname,wigname)
%
% targetlist points to tab-delimited file with one line per target region:  <targname/genename>   <chr>   <start>   <end>
%
% outname = cell array of absolute paths to the output FWB files to be written
%
% wigname = cell array of absolute paths to the wiggle files.
%                      some/all members of the cell arrays can be cell arrays of wiggle files to "OR" together.
%
% Mike Lawrence 2010-09-12

width = 1;

demand_file(targetlist);
targcols = '2 3 4';

java_classpath = '/xchip/cga2/lawrence/cga/trunk/analysis_pipeline/tools/classes';
java_class = 'org.broadinstitute.cga.tools.seq.WigToFwb';

if ~iscell(outname) || ~iscell(wigname)
  error('outname and wigname should both be cell-arrays');
end

if length(outname) ~= length(wigname), error('length(outname)~=length(wigname)'); end
ns = length(outname);
infile = cell(ns,1);
outfile = cell(ns,1);
for s=1:ns
  outfile{s} = outname{s};
  if ischar(wigname{s})
    infile{s} = {wigname{s}};
  elseif iscell(wigname{s})
    infile{s} = wigname{s};
  end
end

% SUBMIT JOBS

did_nothing = true;
all_done = false;
while(~all_done)
  jobs=[];
  cmds = {}; banners = {};
  all_done = true;
  for s=1:ns
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
    banners{end+1,1} = [num2str(s) 'WFB'];
    cmds{end+1,1} = ['-R "rusage[mem=7]" "java -Xmx5g -classpath ' java_classpath ' ' java_class ...
        ' ' targetlist ' ' targcols ' ' outfile{s} ' ' num2str(width)];
    for i=1:length(infile{s}), cmds{end} = [cmds{end} ' ' infile{s}{i}]; end
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
  bwait(jobs);
end

if did_nothing, fprintf('All files are already up-to-date\n');
else fprintf('All files are now up-to-date\n');
end

