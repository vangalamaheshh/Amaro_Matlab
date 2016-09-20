function M = run_MutSig1x(varargin)
% two ways to run:
%   run_MutSig1x(libdir,isetname,maf,outdir)
%   run_MutSig1x(maf,outdir,P)

%%% figure out if "libdir" has been supplied as first argument
if exist(varargin{1},'dir')
  if nargin~=5, error('usage: run_MutSig1.x(libdir,isetname,maf,outdir,P)'); end
  libdir = varargin{1};
  isetname = varargin{2};
  maf = varargin{3};
  outdir = varargin{4};
  P = varargin{5};
elseif exist(varargin{1},'file')
  if nargin<2 || nargin>3, error('usage: run_MutSig1.x(maf,outdir,P)'); end
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
P=impose_default_value(P,'cosmic_file','/xchip/cga1/annotation/db/cosmic/CosmicMutantExport_v54_120711.tsv');
P=impose_default_value(P,'covfile','');

% PREPARE OUTPUT DIRECTORY

if ~exist('isetname','var'), isetname = 'mutsig'; end
outstem = [outdir '/' isetname];  % for classic pipeline, don't add the "."
ede(outdir);

% LOAD DATA

if strcmp(P.covfile,'')
  fprintf('Using generic coverage file\n');
  M = new_load_mutdata(maf,P);
else
  fprintf('Using specified coverage file\n');
  P.mutfile = maf;
  M = load_all_mutation_data2(P);
end

% ANALYSIS 

%try
  P.quick = true;
  M = analyze_mutation_data(M,outstem,P);
%catch me
%  fprintf('ERROR in analyze_mutation_data:\n');
%  disp(me);
%  disp(me.message);
%end

fprintf('Finished, wrote results to %s\n',outdir);

if exist('libdir','var')     % if running in Firehose-deployed mode, then:
  close all                  % close all figures so xvnc can terminate
end

