function M = run_MutSigCV(varargin)
% two ways to run:
%   run_MutSigCV(libdir,isetname,maf,outdir)
%   run_MutSigCV(maf,outdir,P)

%%% figure out if "libdir" has been supplied as first argument
if exist(varargin{1},'dir')
  if nargin~=5, error('usage: run_MutSigCV(libdir,isetname,maf,outdir,P)'); end
  libdir = varargin{1};
  isetname = varargin{2};
  maf = varargin{3};
  outdir = varargin{4};
  P = varargin{5};
elseif exist(varargin{1},'file')
  if nargin<2 || nargin>3, error('usage: run_MutSigCV(maf,outdir,P)'); end
  maf = varargin{1};
  outdir = varargin{2};
  if nargin>=3
    P = varargin{3};
  end
else
  error('NOT FOUND: %s',varargin{1});
end

if exist('libdir','var')
  % set java classpath to make accessible all jars in the libdir
  javaclasspath([libdir;direc([libdir '/*.jar'])]);
end

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'apply_mutation_blacklist',true);
P = impose_default_value(P,'canned_bagels','/xchip/cga/reference/mutsig_params/bagels.v2.mat');
P = impose_default_value(P,'bagelfile',P.canned_bagels);
P = impose_default_value(P,'cosmic_file','/xchip/cga/reference/annotation/db/cosmic/CosmicMutantExport_v54_120711.tsv');
P = impose_default_value(P,'covfile','');
P = impose_default_value(P,'generate_coverage_plot_and_bargraphs',true);
P = impose_default_value(P,'use_corrected_expansion_method',true);

% PREPARE OUTPUT DIRECTORY
if exist('isetname','var')
  outstem = [outdir '/' isetname '.'];
else
  outstem = [outdir '/'];
end
ede(outdir);

% load data
if strcmp(P.covfile,'')
  fprintf('Using generic coverage file\n');
  M = new_load_mutdata(maf,P);
else
  fprintf('Using specified coverage file\n');
  P.mutfile = maf;
  M = load_all_mutation_data2(P); 
  if ~isfield(M,'n_flank')
    M.n_flank = M.n_nonsilent;
    M.n_flank(:) = 0;
  end
  if ~isfield(M,'N_flank_cov')
    M.N_flank_cov = M.N_non_cov;
    M.N_flank_cov(:) = 0;
  end
end

% output version report textfile
if ~P.use_corrected_expansion_method
  ver = 'v0.6';
else
  ver = 'v0.9';
end
save_textfile(['MutSigCV ' ver],[outstem 'MutSig_version.txt']);

% add nsite, npat
M = count_sites_and_patients(M);

% significance analysis
M = estimate_mutrates_CV(M,P);   % (uses P.use_corrected_expansion_method)
P.gene_names = M.gene.name; P.patient_names = M.patient.name;
P.sig_calculation_method = 'projection';
P.null_boost_factor = 1000;
P.projection_bin_cap = 0;
nc = M.TOT-1;
[p score aux] = calculate_significance(M.N_non_cov(:,1:nc,:),M.n_nonsilent(:,1:nc,:),...
                                       {M.mutrate_analysis.X,M.mutrate_analysis.x},P);
% add FDR, sort
g = keep_fields(M.gene,{'name','Nnon','Nsil','Nflank','nnon','npat','nsite','nsil','nflank','nnei','fMLE'});
g.p=p; g.score = score; g.time = aux{1};
g.q=calc_fdr_value(g.p);
g = sort_struct(g,'p');

% write report
fname = [outstem 'sig_genes.txt'];
M.g = rename_field(g,'name','gene');
save_struct(M.g,fname);

fname = [outstem 'results.mat'];
gene = M.gene; g = M.g;
save(fname,'gene','g');

pr(M.g,1:30)

% after the main analysis, do the secondary tasks
%(because these take a long time, and the graphics can crash matlab if X-win not running)

% write maf
save_struct(M.mut,[outstem 'final_analysis_set.maf']);

% summary analyses
if P.generate_coverage_plot_and_bargraphs
  generate_coverage_plot_and_bargraphs(M,outstem,P)
end
save_struct(M.categ,[outstem 'mutcategs.txt']);

%fname = [outstem 'qq_plot.png'];
%print_to_file(fname);  % NOTE, this hangs the terminal indefinitely if XVNC is not connected to the session

fprintf('Finished, wrote results to %s\n',outdir);

if exist('libdir','var')     % if running in Firehose-deployed mode, then:
  close all                  % close all figures so xvnc can terminate
end

