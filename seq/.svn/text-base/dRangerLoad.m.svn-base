function [T N] = dRangerLoad(sample,P)

if ~exist('sample','var'), error('<sample> is required'); end

if ~exist('P','var'), P=[]; end
if ~isempty(P) && ~isstruct(P), error('Second parameter should be a P struct'); end

P = impose_default_value(P,'cancer_sample',true);
P = impose_default_value(P,'dRangerPreprocess_output_dir_suffix','dR');
P = impose_default_value(P,'dRangerPrepare_output_file','all.weird.joined.mat');

% eventually will be:
%P = impose_default_value(P,'dRangerPreprocess_output_dir_suffix','dR2');
%P = impose_default_value(P,'dRangerPrepare_output_file','stringent_pairs.mat');

f = [P.dRangerPreprocess_output_dir_suffix '/' P.dRangerPrepare_output_file];

if P.cancer_sample
  T = load(['/xchip/tcga_scratch/lawrence/' sample '/tumor_' f]);
  T = T.x;
  N = load(['/xchip/tcga_scratch/lawrence/' sample '/normal_' f]);
  N = N.x;
else
  T = load(['/xchip/tcga_scratch/lawrence/' sample '/sample_' f]);
  T = T.x;
  N = [];
end
