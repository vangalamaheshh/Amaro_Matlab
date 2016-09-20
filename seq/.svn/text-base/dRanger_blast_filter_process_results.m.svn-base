function X = dRanger_blast_filter_process_results(sample,P,X);
% dRanger_blast_filter_process_results(sample,P,[,X]);
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

try

fprintf('\ndRanger_blast_filter_process_results\n\tsample = %s\n',sample);

% directories and filenames

short_sample = sample_to_short_sample(sample);
direc = ['/xchip/tcga_scratch/lawrence/' sample];
outdir = direc;
tempdir = [direc '/' P.results_name '_' P.blast_filter_working_directory_suffix];
if ~exist(tempdir,'dir'), error('Cannot find directory %s',tempdir); end
prefix = [tempdir '/blast_' ];
intable = [prefix 'intable.txt'];
outtable = [prefix 'outtable.txt'];
hitstem = [prefix 'hit_batch_'];
fname = [direc '/' P.results_name '.txt'];

% load initial dRanger results

if ~exist('X','var')
  fprintf('Loading dRanger results\n');
  if ~exist(fname,'file'), error('Can''t find file %s',fname);end
  X = load_struct(fname);
end

X = make_numeric(X,{'chr1','chr2','min1','min2','max1','max2','normreads'});

if isfield(X,'filterB')
   fprintf('\n%s already has a filterB column.\n');
   fprintf('re-running anyways\n');
end

if P.blast_somatic_only
  X1 = reorder_struct(X,X.normreads==0);
else
  X1 = X;
end

nx1 = slength(X1);

% check to make sure all blast results files are present,
% and find out the date of the newest one

nf = ceil(nx1/P.max_queries_per_blast_batch);

all_present = true;
latest_datenum = [];
for f=1:nf
  hfname = [tempdir '/blast_hit_batch_' num2str(f)];
  d = dir(hfname);
  if isempty(d), fprintf('Not found: %s\n',hfname); all_present = false; continue; end
  if isempty(latest_datenum) | d.datenum>latest_datenum, latest_datenum = d.datenum; end
end
if ~all_present
  error('Aborting: not all required files are present.  Re-run dRanger_blast_filter.');
end

% check for already-existing outtable

d = dir(outtable);
if ~isempty(d)
  fprintf('\n%s already exists.\n', outtable);
  if d.datenum < latest_datenum
    fprintf('  But it is older than at least one of the blast results files.\n');
    fprintf('  Deleting and refreshing.\n');
    delete(outtable);
  else
    fprintf('  And it doesn''t need to be refreshed.\n');
  end
end

% process_results script

scriptname = '/home/radon00/lawrence/CancerGenomeAnalysis/trunk/matlab/seq/dRanger_blast_filter_process_results.pl';

% dRanger_blast_filter_process_results.pl
% input:
%    --> filestem
%    --> table of rearrangements
%           jump# file# chr1 mn1 mx1 chr2 mn2 mx2
% otuput:
%    --> table of rearrangements
%           jump# hits1 self1 other1 hits2 self2 other2
%                 (hits = total # of hits)
%                 (self = number of hits to self)
%                 (other = number of hits to pairmate)


pos1 = round((X1.min1+X1.max1)/2);
pos2 = round((X1.min2+X1.max2)/2);
J = []; J.jump = (1:nx1)';
J.file = ceil(J.jump/P.max_queries_per_blast_batch);
J.chr1 = X1.chr1; J.mn1 = pos1-P.blast_radius; J.mx1 = pos1+P.blast_radius;
J.chr2 = X1.chr2; J.mn2 = pos2-P.blast_radius; J.mx2 = pos2+P.blast_radius;

save_struct(J,intable,'no_headers');
banner = [short_sample 'DBFPP'];
cmd = ['"perl ' scriptname ' ' intable ' ' hitstem ' ' outtable '"'];
n_tries = 0;
while(1)
  if exist(outtable,'file'), break; end
  n_tries = n_tries + 1;
  bwait(bsub(cmd,banner));
  if exist(outtable,'file'), break; end
  if n_tries == 10, error('DBFPP perlscript failed 10 times!\n'); end
  fprintf('DBFPP perlscript failed... trying again\n');
end

fprintf('Retrieving processed results from %s\n', outtable);
K = load_struct(outtable,'%*f%f%f%f%f%f%f',0);
K = rename_field(K,{'col1','col2','col3','col4','col5','col6'},...
  {'hits1','self1','other1','hits2','self2','other2'});

z = nan(slength(X),1);
X.hits1=z; X.self1=z; X.other1=z;
X.hits2=z; X.self2=z; X.other2=z;

if P.blast_somatic_only
  idx = find(X.normreads==0);
else
  idx = 1:slength(X);
end

X.hits1(idx) = K.hits1; X.self1(idx) = K.self1; X.other1(idx) = K.other1;
X.hits2(idx) = K.hits2; X.self2(idx) = K.self2; X.other2(idx) = K.other2;

X.filterB = 1*(X.other1>0 | X.other2>0);

if P.blast_somatic_only
  X.filterB(X.normreads>0) = -1;
end

% SAVE RESULTS
fname = [direc '/' P.results_name '.txt'];
fprintf('Saving results to %s\n',fname);
save_struct(X,fname);

fprintf('Done!\n');

catch me, excuse(me); end
