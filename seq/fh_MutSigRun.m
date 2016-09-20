function M = fh_MutSigRun(libdir, varargin)
% fh_MutSigRun
%
% performs analysis and report generation for MutSig
%
% requires the following parameters:
%    -i <individual_set_id>
%    -b <build> = 'hg18', 'hg19', etc.
%    -bd <build_dir> = genome reference sequence directory with chr*.txt files and build_info.txt
%    -maf <maffile> = maf file produced by fh_MutSigPreprocess
%    -cov <covfile> = coverage.mat file produced by fh_MutSigCoverageGather
%    -pat <patlist> = patients.txt file produced by fh_MutSigPreprocess
%    -ppr <preprocessing_report> = mutation_preprocessing_report.txt file produced by fh_MutSigPreprocess
%    -cat <categories_list> = filename of list of final categories for MutSig (e.g. CpG,other CG,AT,Indel,Null)
%    -genes <genelist> = list of genes to consider during MutSig
%    -genesets <genesetfile> = list of genesets to consider in genesets analysis
%    -cosmic <cosmicfile> = cosmic file to use in cosmic analysis
%    -refseq <refseqfile> = refseq mat file to use
%    -jobcount <jobcount> = how many scatter-jobs have been dispatched?
%    -whichjob <whichjob> = which job am I?  1-<jobcount> for scatter, or 0 for gather
%
% optional parameters
%    -p <paramfile> (optional) = file with column1=key and column2=value
%
% outputs many files, including report

% Mike Lawrence 2010-09-28

fprintf('fh_MutSigRun\n');

if ~exist('libdir','var') || ~ischar(libdir) || ~exist(libdir,'dir')
  error('First parameter needs to be libdir path');
end

fprintf('libdir = %s\n',libdir);
fprintf('%s\n',varargin{:});
fprintf('\n');

% set java classpath to make accessible all jars in the libdir
javaclasspath([libdir;direc([libdir '/*.jar'])]);

% FIREHOSE PARAMETERS
required_flds = {'b','bd','maf','cov','genesets','cosmic','refseq'};
optional_flds = {'i','pat','genes','cat','contextdir','p','ppr','jobcount','whichjob'};
args = handle_args([required_flds,optional_flds],varargin);
require_fields(args,required_flds);

% INDIVIDUAL_SET_ID
if isempty(args.i)
  fprintf('No individual_set_id provided: will use "MutSig" as outputfile prefix.\n');
  args.i = 'MutSig';
end

% MUTSIG PARAMETERS
P = [];
P = process_params_file(P,args.p);

P = impose_default_value(P,'patlist',args.pat);
P = impose_default_value(P,'covfile',args.cov);
P = impose_default_value(P,'mutfile',args.maf);
P = impose_default_value(P,'catfile',args.cat);
P = impose_default_value(P,'genelist',args.genes);
P = impose_default_value(P,'build',args.b);
P = impose_default_value(P,'contextdir',args.contextdir);
P = impose_default_value(P,'geneset_collection_file',args.genesets);
P = impose_default_value(P,'build_cosmic',args.b);
P = impose_default_value(P,'cosmic_file',args.cosmic);
P = impose_default_value(P,'build_refseq',args.refseq);
P = impose_default_value(P,'build_dir',args.bd);
P = impose_default_value(P,'build_genome_region',args.bd);
P = impose_default_value(P,'quick',false);
P = impose_default_value(P,'output_per_exon_coverage',true);
P = impose_default_value(P,'output_per_exon_mutations',false);
P = impose_default_value(P,'output_final_mutation_list',true);
P = impose_default_value(P,'output_per_gene_coverage',true);
P = impose_default_value(P,'output_per_gene_mutation_counts',true);
P = impose_default_value(P,'output_sample_sig_gene_table',true);
P = impose_default_value(P,'force_recalculate_maf_simplename_fields',true);
P = impose_default_value(P,'mutation_preprocessing_report_file',args.ppr);
P = impose_default_value(P,'mutation_filtering_report_file',['./' args.i '.mutation_filtering_report.txt']);
P = impose_default_value(P,'print_report', true);

% CHECK PARALLELIZATION
if ~isempty(args.jobcount) && ~isempty(args.whichjob)
  jobcount = str2double(args.jobcount);
  whichjob = str2double(args.whichjob);
  if isnan(jobcount) || isnan(whichjob) || ~isround(jobcount) || ~isround(whichjob) ||...
    jobcount<1 || whichjob<0 || whichjob>jobcount
    error('Improper values: <jobcount>=%s, <whichjob>=%s',args.jobcount,args.whichjob);
  end
  P.using_scatter_gather = true;
  P.scatter_gather_jobcount = jobcount;
  P.scatter_gather_whichjob = whichjob;
  P.is_scatter = (whichjob>0);
  P.is_gather = (whichjob==0);
elseif ~isempty(args.jobcount) || ~isempty(args.whichjob)
  error('<jobcount> and <whichjob> cannot be specified except together');
else
  P.using_scatter_gather = false;
end

disp(P)

% initialize ReferenceInfoObj from the genome_dir
try
  ReferenceInfoObj.init(P.build_dir);
catch me
  disp(me); disp(me.message);
  fprintf('WARNING: Failed to initialize ReferenceInfoObj from P.build_dir = %s\n',args.bd);
  if isempty(args.bd), fprintf('\tbecause P.build_dir was unspecified!\n'); end
end

% output directory
outdir = '.';

% MUTSIG
M = mutsig_analysis_and_report(args.i,outdir,P);

close all    % close all figures so xvnc can terminate

