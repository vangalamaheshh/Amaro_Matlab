function S = load_S_perlane_data(sample,P)
% S = load_S_perlane_data(sample,P)
%
% Mike Lawrence 2009-08-04

if ~exist('sample','var'), error('<sample> is required'); end
if iscell(sample), error('Does not support multiple samples.'); end

if ~exist('P','var'), P=[]; end
if ~isempty(P) && ~isstruct(P), error('Second parameter should be a P struct'); end

P = impose_default_value(P,'InsertSizeByLane_output_dir_suffix','isz');

basedir = '/xchip/tcga_scratch/lawrence/';
S = process_isz_data([basedir sample '/sample_' P.InsertSizeByLane_output_dir_suffix]);

S.lanes = load_struct([basedir sample '/sample.bam.lanelist'],'%f%s',0);
fname = [basedir sample '/sample.bam.lanetable'];
if exist(fname,'file')
  tmp = load_struct_specify_numeric_cols(fname,[1]);
  S.lanes = merge_structs({S.lanes,tmp});
end

% colors
S.color = zeros(S.nlanes,3); % black

