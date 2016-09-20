function M = run_MutSigS2N(varargin)
% two ways to run:
%   run_MutSigS2N(libdir,isetname,maf,outdir)
%   run_MutSigS2N(maf,outdir,P)

%%% figure out if "libdir" has been supplied as first argument
if exist(varargin{1},'dir')
  if nargin~=5, error('usage: run_MutSigS2N(libdir,isetname,maf,outdir,P)'); end
  libdir = varargin{1};
  isetname = varargin{2};
  maf = varargin{3};
  outdir = varargin{4};
  P = varargin{5};
elseif exist(varargin{1},'file')
  if nargin<2 || nargin>3, error('usage: run_MutSigS2N(maf,outdir,P)'); end
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
P=impose_default_value(P,'cosmic_file','/xchip/cga/reference/annotation/db/cosmic/CosmicMutantExport_v54_120711.tsv');
P=impose_default_value(P,'covfile','');

% prepare output directory
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
end

% output version report textfile
save_textfile('MutSigS2N v0.3',[outstem 'MutSig_version.txt']);

% analysis
M = simple_calc_rates_v2(M,P);
pr(M.g,1:30)

fname = [outstem 'sig_genes.txt'];
M.g = rename_field(M.g,'name','gene');
save_struct(M.g,fname);

fname = [outstem 'mutcategs.txt'];
save_struct(M.categ,fname);

fname = [outstem 'results.mat'];
gene = M.gene; g = M.g;% rate=keep_fields(M,grep('_c$',fieldnames(M)));
save(fname,'gene','g');

fname = [outstem 'qq_plot.png'];
%print_to_file(fname);  % NOTE, this hangs the terminal indefinitely if XVNC is not connected to the session

fprintf('Finished, wrote results to %s\n',outdir);

if exist('libdir','var')     % if running in Firehose-deployed mode, then:
  close all                  % close all figures so xvnc can terminate
end

