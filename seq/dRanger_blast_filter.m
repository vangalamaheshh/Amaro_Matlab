function X = dRanger_blast_filter(sample,P,X);
% dRanger_blast_filter(sample,P[,X])
%
% Mike Lawrence 2009

if ~exist('sample','var'), error('<sample> is required'); end

if ~exist('P','var'), P=[]; end
if ~isempty(P) && ~isstruct(P), error('Second parameter should be a P struct'); end

P=impose_default_value(P,'blast_radius',250);
P=impose_default_value(P,'blast_somatic_only',true);
P=impose_default_value(P,'blast_filter_working_directory_suffix','dbf2');
P=impose_default_value(P,'max_queries_per_blast_batch',20);
P=impose_default_value(P,'results_name','dRanger_results');

batch_blast_script = '/home/radon01/lawrence/CancerGenomeAnalysis/trunk/matlab/seq/batch_blast.csh';

try

direc = ['/xchip/tcga_scratch/lawrence/' sample];
outdir = direc;

fprintf('\ndRanger_blast_filter\n\tsample = %s\n',sample);
short_sample = sample_to_short_sample(sample);

% load initial dRanger results

fname = [direc '/' P.results_name '.txt'];

if ~exist('X','var')
  fprintf('Loading dRanger results\n');
  if ~exist(fname,'file'), error('Can''t find file %s',fname);end
  X = load_struct(fname);
end

X = make_numeric(X,{'chr1','chr2','min1','min2','max1','max2'});

if isfield(X,'filterB')
  fprintf('%s already has a filterB column.\n');
  fprintf('re-running anyways\n');
end

nx = slength(X);

% temporary directory

tempdir = [direc '/' P.results_name '_' P.blast_filter_working_directory_suffix];
if ~exist(tempdir,'dir'), mkdir(tempdir), end

% make list of blast queries

did_nothing = true;
fprintf('Generating blast queries... ');


if P.blast_somatic_only
  X = make_numeric(X,'normreads');
  X1 = reorder_struct(X,X.normreads==0);
else
  X1 = X;
end

nx1 = slength(X1);

query_files={}; nf=0;
q=''; nq=0;
for i=1:nx1
  if ~mod(i,100), fprintf('%d/%d ',i,nx1); end
  chr1 = X1.chr1(i); pos1 = round((X1.min1(i)+X1.max1(i))/2);
  chr2 = X1.chr2(i); pos2 = round((X1.min2(i)+X1.max2(i))/2);
  mn1 = pos1-P.blast_radius; mx1 = pos1+P.blast_radius; probe1 = genome_region(chr1,mn1,mx1);
  mn2 = pos2-P.blast_radius; mx2 = pos2+P.blast_radius; probe2 = genome_region(chr2,mn2,mx2);
  q = [q '>jump' num2str(i) '.end1' char(10) probe1 char(10)...
         '>jump' num2str(i) '.end2' char(10) probe2 char(10)];
  nq = nq + 1;
  if nq==P.max_queries_per_blast_batch || i==nx1
    nf = nf + 1;
    fname = [tempdir '/blast_query_batch_' num2str(nf)];
    if ~exist(fname,'file')
      save_textfile(q,fname);
      did_nothing = false;
    end
    query_files = [query_files; fname];
    q=''; nq=0;
  end
end
if did_nothing
  fprintf('All query files already existed\n');
end

% farm out blast batches

hit_files = regexprep(query_files,'query','hit');

fprintf('\nSubmitting blast batches to LSF\n');
did_nothing = true;
all_done = false;
while(~all_done)
  all_done = true;
  jobs=[];
  for f=1:nf
    qfile = query_files{f};
    hfile = hit_files{f};
    dhfile = dir(hfile);
    if ~isempty(dhfile) & dhfile.bytes > 200, continue; end
    all_done = false;
    did_nothing = false;
    banner = [short_sample 'BLST' num2str(f)];

    cmd = ['"' batch_blast_script ' ' qfile ' ' hfile '"'];

%%%%%%%%%%%%%%% original approach:
%    cmd = ['"blastall -p blastn -d /xchip/tcga/gbm/analysis/lawrence/genome/hg18/orig/hg18.fa ' ...
%      '-m8 -e1e-50 -nT -FF < ' qfile ' > ' hfile '_partial;mv ' hfile '_partial ' hfile '"'];
%%%%% this approach had a bug: the "mv" command would get executed even if the blastall got aborted.

    jobs = [jobs; bsub(cmd,banner)];
  end
  if ~all_done
    % wait for jobs to finish
    fprintf('\nWaiting for blast batches to finish\n');
    bwait(jobs);
  end
end % while(~all_done)

if did_nothing
  fprintf('All results files already existed\n');
end

X = dRanger_blast_filter_process_results(sample,P,X);

catch me, excuse(me); end
