function F = crawl_firehose(workspace_id,outdir)
% Mike Lawrence 2011-08-24

if ~exist('outdir','var'), outdir = '/tmp/riker'; end
ensure_dir_exists(outdir);

surveyfhscript = '/xchip/cga2/lawrence/cga/trunk/matlab/seq/surveyfh.pl';
outfile = [outdir '/' workspace_id '.survey.txt'];

result = system(['perl ' surveyfhscript ' ' workspace_id ' ' outfile]);
F = load_struct_noheader(outfile,{'individual_set_id','sample_id','exp','bampath'});
F = parse_in(F,'sample_id','^(.*)-(Tumor|Normal)$',{'indiv','tn'});
F.tn = lower(F.tn);

%% returns F with
%      individual_set_id
%      sample_id
%   *  indiv = individual name
%   *  exp = capture or wgs
%   *  tn = tumor or normal
%   *  bampath = picard BAM path

