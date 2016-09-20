function [T N] = load_GC_perlane_data(sample,P)
% [T N] = load_GC_perlane_data(sample,P)
%
% Mike Lawrence 2009-07-30

if ~exist('sample','var'), error('<sample> is required'); end
if iscell(sample), error('Does not support multiple samples.'); end

if ~exist('P','var'), P=[]; end
if ~isempty(P) && ~isstruct(P), error('Second parameter should be a P struct'); end

P = impose_default_value(P,'GCContentByLane_output_dir_suffix','gc');

basedir = '/xchip/tcga_scratch/lawrence/';
T = process_gc_data([basedir sample '/tumor_' P.GCContentByLane_output_dir_suffix]);
N = process_gc_data([basedir sample '/normal_' P.GCContentByLane_output_dir_suffix]);

T.lanes = load_struct([basedir sample '/tumor.bam.lanelist'],'%f%s',0);
fname = [basedir sample '/tumor.bam.lanetable'];
if exist(fname,'file')
  tmp = load_lanetable(fname);
  T.lanes = merge_structs({T.lanes,tmp});
end

N.lanes = load_struct([basedir sample '/normal.bam.lanelist'],'%f%s',0);
fname = [basedir sample '/normal.bam.lanetable'];
if exist(fname,'file')
  tmp = load_lanetable(fname);
  N.lanes = merge_structs({N.lanes,tmp});
end


