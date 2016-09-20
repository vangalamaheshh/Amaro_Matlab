function fh_MutLego(libdir, varargin)
% fh_MutLego
%
% makes lego plots  and report for individual_set based on maf file
%
% requires the following parameters:
%    -i <individual_set_id>
%    -maf <maffile> = maf file produced by fh_MutSigPreprocess
%    -cov <covfile> = coverage.mat file produced by fh_MutSigCoverageGather
%
% optional parameters
%    -p <paramfile> (optional) = file with column1=key and column2=value
%
% outputs many files, including report html

% Chip Stewart 2012-08-14

fprintf('fh_MutLego\n');
fprintf('libdir = %s\n',libdir);
fprintf('%s\n',varargin{:});
fprintf('\n');

% FIREHOSE PARAMETERS
required_flds = {'maf'};
optional_flds = {'i','cov'};
args = handle_args([required_flds,optional_flds],varargin);
require_fields(args,required_flds);

% MUTSIG PARAMETERS
P = [];
P = process_params_file(P,args.p);

P = impose_default_value(P,'covfile',args.cov);
P = impose_default_value(P,'mutfile',args.maf);
P = impose_default_value(P,'print_report', true);

disp(P)

% individual_set_id and output directory

if isempty(args.i)
  fprintf('No individual_set_id provided: will use "lego" as outputfile prefix.\n');
  args.i = 'lego';
end

outdir = '.';

% MUTSIG
plotMutationSpectrumCategLegos(P.mutfile,P.outdir,P);

close all    % close all figures so xvnc can terminate

