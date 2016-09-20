function X = dRanger_blast_filter_process_results2(jobid,X,outdir,tempdir,max_queries_per_blast_batch,blast_radius,max_hits)

if ~exist('jobid','var'), error('jobid is required'); end
if ~exist('X','var'), error('X is required'); end
if ~exist('outdir','var'), error('outdir is required'); end

if ~exist('max_queries_per_blast_batch','var'), max_queries_per_blast_batch = 25; end
if ~exist('blast_radius','var'), blast_radius = 250; end
if ~exist('max_hits','var'), max_hits = 10; end
if ~exist('tempdir','var'), tempdir = '/xchip/tcga_scratch/lawrence/tmp'; end

fprintf('dRanger_blast_filter_process_results2\n\toutdir = %s\n\tjobid = %s\n',outdir,jobid);

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

pos1 = round((X.min1+X.max1)/2);
pos2 = round((X.min2+X.max2)/2);
nx = slength(X);
J = []; J.jump = (1:nx)';
J.file = ceil(J.jump/max_queries_per_blast_batch);
J.chr1 = X.chr1; J.mn1 = pos1-blast_radius; J.mx1 = pos1+blast_radius;
J.chr2 = X.chr2; J.mn2 = pos2-blast_radius; J.mx2 = pos2+blast_radius;

prefix = [tempdir '/blast_' jobid '_'];
intable = [prefix 'intable.txt'];
outtable = [prefix 'outtable.txt'];
hitstem = [prefix 'hit_batch_'];
save_struct(J,intable,'no_headers');
cmd = ['perl ' scriptname ' ' intable ' ' hitstem ' ' outtable];
job = bsub(['"' cmd '"']);
bwait(job);

K = load_struct(outtable,'%f%f%f%f%f%f%f',0);
K = rename_field(K,{'col1','col2','col3','col4','col5','col6','col7'},...
  {'jump','hits1','self1','other1','hits2','self2','other2'});

X = merge_structs({X,K});
X.filterB = zeros(nx,1);
X.filterB(X.hits1>max_hits | X.hits2>max_hits) = 1;
X.filterB(X.other1>0 | X.other2>0) = 2;
X.filterB(X.self1<1 | X.self2<1) = NaN;   % BLAST failed
