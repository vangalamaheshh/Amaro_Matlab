function SegSeqRun(sample,P)
% SegSeqRun(sample)
%
% sample: e.g. ov/0751/wgs
%
% Mike Lawrence 2009-08-03
% based on code from Gordon Saksena

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'use_LSF',true);

if ~exist('sample','var'), error('<sample> required'); end
if iscell(sample), error('Multiple samples not supported'); end
if ~contains(sample,'wgs') || contains(sample,'-') || ~contains(sample,'/')
  error('<sample> should be of form ov/0751/wgs');
end

name2 = upper(regexprep(sample,'/','-'));
name2 = regexprep(name2,'-WGS$','');

% make sure preprocessing directories are linked in from "ng" dir

bd = '/xchip/tcga_scratch/lawrence';
ng = '/xchip/tcga_scratch/ng';
cn = [ng '/' name2 '/wgs/cn/'];
if ~exist(cn,'dir'), mkdir(cn); end
force_link([bd '/' sample '/tumor_SS'],cn);
force_link([bd '/' sample '/normal_SS'],cn);

% run SegSeq
if P.use_LSF
  cmd = ['-q hugemem -R "rusage[matlab=1]" /home/radon00/lawrence/CancerGenomeAnalysis/trunk/matlab/seq/SegSeqRun.csh ' name2];
  banner = [name2 '_SS'];
  bwait(bsub(cmd,banner));
else
  cmd = ['/home/radon00/lawrence/CancerGenomeAnalysis/trunk/matlab/seq/SegSeqRun.csh ' name2];
  system(cmd);
end

% link results into "lawrence" directory

SegSeq_link_best(sample);

