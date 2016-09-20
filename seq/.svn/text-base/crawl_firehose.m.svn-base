function F = crawl_firehose(fdir)
% Mike Lawrence 2009-09-15

if ~exist('fdir','var'), fdir = '/xchip/cga1/firehose_output/trunk/Individual'; end
jdir = '/xchip/cga/gdac-prod/genepattern/jobResults';

wd = pwd;
cd(fdir);

flinks_t = split(ls('-ltr','*/*/bam/raw/tumor.bam'),char(10));
flinks_n = split(ls('-ltr','*/*/bam/raw/normal.bam'),char(10));
flinks = [flinks_t;flinks_n];
flinks = remove_ansi(flinks);
cd(wd);

jcrawl_script = '/home/radon00/lawrence/CancerGenomeAnalysis/trunk/matlab/seq/jcrawl.csh';
[tmp jj] = system([jcrawl_script ' ' jdir]);
jlinks = split(jj,char(10));
jlinks = remove_ansi(jlinks);

F = parse(flinks,[' (\S*/.*.bam).*(' jdir '/.*.bam)'], {'ifile','jfile'});
F.ifile = regexprep(F.ifile,'^(.*)$',[fdir '/$1']);

J = parse(jlinks,[' (\S*/.*.bam) -> (/.*.bam)'], {'jfile','pfile'});
J.jfile = regexprep(J.jfile,'^(.*)$',[jdir '/$1']);

tmp = parse(F.ifile,'Individual/([^/]*)/([^/]*)/bam/raw/([^\.]*)\.bam', {'indiv','exp','tn'});
F = merge_structs({F,tmp});

F.bampath = map_across(F.jfile,J.jfile,J.pfile);

%% returns F with
%      (ifile)
%      (jfile)
%   *  indiv = individual name
%      exp = capture or wgs
%   *  tn = tumor or normal
%   *  bampath = picard BAM path

