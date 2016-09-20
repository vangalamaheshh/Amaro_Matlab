function fh_MutSigS2NRun(libdir, maf, blacklist, isetname, build)
% fh_run_MutSig_S2N
%
% Performs analysis for MutSig_S2N
%
% Requires the following parameters:
%    <libdir>    Module directory where executables reside.
%    <maf>       Concatenated MAF file path.
%    <blacklist> Path to blacklist file.
%    <build>     Build (e.g. hg18, hg19, etc.).
%
% Outputs significant genes file.

% Dan DiCara/Mike Lawrence 2012-08-28

fprintf('fh_MutSigS2NRun\n');
fprintf('libdir = %s\n',libdir);
fprintf('maf = %s\n',maf);
fprintf('blacklist = %s\n',blacklist);
fprintf('individual set name = %s\n',isetname);
fprintf('build = %s\n',build);
fprintf('\n');

% output version report textfile
MUTSIG_VERSION = 'S2N'
outname = [isetname '.MutSig_version.txt'];
save_textfile(['MutSig v' MUTSIG_VERSION],outname);

% set java classpath to make accessible all jars in the libdir
if ~exist(libdir, 'dir')
    error('libdir cannot be found: %s', libdir)
end

% Check that MAF file exists
if ~exist(maf,'file')
    error('MAF cannot be found: %s', maf)
end

% MutSig S2N parameters
if ~exist(blacklist, 'file')
    error('blacklist cannot be found: %s', blacklist)
end
P = struct();
P.mutation_blacklist = blacklist;
P.build= build;
P = impose_default_value(P,'force_recalculate_maf_simplename_fields',true);

% Output directory - set to current working directory.
outdir = '.';

% Run MutSig S2N
run_MutSig_S2N(libdir, isetname, maf, outdir, P);

 % close all figures so xvnc can terminate
close all
