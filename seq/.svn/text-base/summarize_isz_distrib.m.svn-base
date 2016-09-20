function d = summarize_isz_distrib(sample,tn,P)
% summarize_isz_distrib(sample,tn)
%
% Given the results of InsertSizeByLane,
% summarizes the distribution by collapsing similar lanes (groups = libraries)
% and returning a matrix of the following form:
%
% d = insert-size distribution properties
%          each row:  frac mean std
%          default =  0.5  250  15
%                     0.5  400  20
%
% Mike Lawrence 2009-09-10

if ~exist('sample','var'), error('<sample> is required'); end
if ~exist('tn','var'), error('<tn> is required'); end
if iscell(sample), error('Does not support multiple samples.'); end
if ~exist('P','var'), P=[]; end
if ~isempty(P) && ~isstruct(P), error('Second parameter should be a P struct'); end
P = impose_default_value(P,'InsertSizeByLane_output_dir_suffix','isz');

P = impose_default_value(P,'mode_match_margin_percent',5);
P = impose_default_value(P,'std_match_margin_percent',20);

switch lower(tn(1))
  case 't', tn = 'tumor';
  case 'n', tn = 'normal';
  case 's', tn = 'sample';
  otherwise error('tn must be "tumor" or "normal" or "sample"');
end

basedir = '/xchip/tcga_scratch/lawrence/';
D = process_isz_data([basedir sample '/' tn '_' P.InsertSizeByLane_output_dir_suffix]);

d = [];
idx = D.good;
used = false(length(idx),1);
for j=1:length(idx)
  if (used(j)), continue; end
  i = idx(j);
  k = find(~used & abs(100*(D.mode(i)-D.mode(idx))./D.mode(i))<=P.mode_match_margin_percent & ...
          abs(100*(D.std(i)-D.std(idx))./D.std(i))<=P.std_match_margin_percent);
  d = [d; sum(D.cov(idx(k))) mean(D.mode(idx(k))) mean(D.std(idx(k)))];
  used(k) = true;
end
d(:,1) = d(:,1) / sum(d(:,1));
[tmp mi] = min(d(:,1));
d(mi,1) = 1-sum(d(setdiff(1:size(d,1),mi),1));
