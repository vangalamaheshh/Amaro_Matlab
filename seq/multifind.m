function success = multifind(B,dirname,freeze)
% 
%  multifind
%
%  Input:
%     B = list of regions to pull out sequences from
%         fields: chr (numeric), start, end
%
%     freeze = number of data freeze to use
%
%     dirname = directory for writing *.out files
%
%  Mike Lawrence 2009-02-04

if ~exist('B','var'), error('Must specify regions'); end
if ~exist('dirname','var'), error('Must specify output directory'); end
if ~exist('freeze','var'), freeze=4; end

require_fields(B,{'chr','start','end'});

% WRITE TARGET FILE

targ_file = [dirname '/targets.txt'];
if ~exist(dirname,'dir'), mkdir(dirname), end;
X=[];
X.chr = B.chr;
X.start = B.start;
X.end = B.end;
save_struct(X,targ_file,'no_headers');

% CALL PERLSCRIPT IN BATCH

script = '/home/radon00/lawrence/CancerGenomeAnalysis/trunk/matlab/seq/multifind.pl';

% submit jobs

deliverables={}; jobs={};
for tn={'tumor','normal'}
  input_dir = sprintf('/xchip/tcga/gbm/analysis/lawrence/wgs/freeze%d/transloc/%s',freeze,tn{1});
  d = dir(input_dir);
  for i=1:length(d)
    name = d(i).name;
    if isempty(regexp(name,'^\d\d\d[ABCDE]$')), continue; end
    in_file = [input_dir '/' name];
    out_file = [dirname '/' name '.out'];
    cmd = ['"perl ' script ' ' targ_file ' ' in_file ' ' out_file '"'];
    job_no = bsub(cmd);
    deliverables = [deliverables; out_file];
    jobs = [jobs; job_no];
  end
end

% wait until all jobs finish

while(1)
  j = bjobs();
  if slength(j)==0, break; end
  still_running = reorder_struct(j,find(ismember(j.jobid,jobs)));
  if slength(still_running)==0, break; end
  pause(60);
end

% check whether all output files were created

success = true;
for i=1:length(deliverables)
  if ~exist(deliverables{i},'file'), success = false; end
end

if success, fprintf('Multifind succeeded! =)\n');
else fprintf('Multifind failed =(\n'); end
