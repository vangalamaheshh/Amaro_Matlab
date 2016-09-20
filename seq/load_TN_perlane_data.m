function [T N] = load_TN_perlane_data(samples,P)
% [T N] = load_TN_perlane_data(samples,P)
%
% Mike Lawrence 2009

if ~exist('samples','var'), error('<samples> is required'); end
if ~iscell(samples), samples = {samples}; end


if ~exist('P','var'), P=[]; end
if ~isempty(P) && ~isstruct(P), error('Second parameter should be a P struct'); end

basedir = '/xchip/tcga_scratch/lawrence/';

for i=1:length(samples)
  sample = samples{i};
  fprintf('\nLOADING SAMPLE %s\n', sample);
  T{i} = load_and_process_isz_data([basedir sample '/tumor_isz'],[basedir sample '/tumor.bam.lanetable'],P);
  N{i} = load_and_process_isz_data([basedir sample '/normal_isz'],[basedir sample '/normal.bam.lanetable'],P);
end

if length(samples)==1
  T = T{1};
  N = N{1};
end
