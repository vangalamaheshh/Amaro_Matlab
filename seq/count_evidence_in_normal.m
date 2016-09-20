function [n_ref n_mut n_other1 n_other2 t_ref t_mut t_other1 t_other2]...
  = count_evidence_in_normal(M,samples,first_idx,last_idx,blacklist)
% Mike Lawrence 2009-11-05

if ~exist('blacklist','var'), blacklist = 'none'; end

require_fields(M,{'chr','start','ref_allele','tum_allele1','tum_allele2','normal_barcode','tumor_barcode'});
require_fields(samples,{'normal_barcode','normal_bam','tumor_barcode','tumor_bam'});

M.chrn = convert_chr(M.chr);
M = make_numeric(M,{'start'});

if ~exist('first_idx','var'),first_idx=1;end
if ~exist('last_idx','var'),last_idx=slength(M);end

base('ACGTN')=1:5;

M.normal_bam = map_across(M.normal_barcode,samples.normal_barcode,samples.normal_bam);
idx = find(cellfun('isempty',M.normal_bam));
if ~isempty(idx)
  fprintf('Don''t know how to find BAMs for the following normals:\n');
  disp(unique(M.normal_barcode(idx)));
  error('Please add them to samples.');
end
M.tumor_bam = map_across(M.tumor_barcode,samples.tumor_barcode,samples.tumor_bam);
idx = find(cellfun('isempty',M.tumor_bam));
if ~isempty(idx)
  fprintf('Don''t know how to find BAMs for the following tumors:\n');
  disp(unique(M.tumor_barcode(idx)));
  error('Please add them to samples.');
end

try


n_ref = nan(slength(M),1);
n_mut = nan(slength(M),1);
n_other1 = nan(slength(M),1);
n_other2 = nan(slength(M),1);
t_ref = nan(slength(M),1);
t_mut = nan(slength(M),1);
t_other1 = nan(slength(M),1);
t_other2 = nan(slength(M),1);
for i=first_idx:last_idx, fprintf('%d/%d ',i,slength(M));
  if strcmp(M.ref_allele{i},M.tum_allele1{i}), newbase = M.tum_allele2{i}; else newbase = M.tum_allele1{i}; end
  [R B S] = pull_from_bam(M.normal_bam{i},M.chrn(i),M.start(i),M.start(i),struct('quiet',1,'blacklist',blacklist));
  bidx = find(B(:,4)==M.start(i));
  n_ref(i) = sum(B(bidx,1)==base(S));
  n_mut(i) = sum(B(bidx,1)==base(newbase));
  others = setdiff(1:4,base([newbase S]));
  n_other1(i) = sum(B(bidx,1)==others(1));
  n_other2(i) = sum(B(bidx,1)==others(2));
  [R B S] = pull_from_bam(M.tumor_bam{i},M.chrn(i),M.start(i),M.start(i),struct('quiet',1,'blacklist',blacklist));
  bidx = find(B(:,4)==M.start(i));
  t_ref(i) = sum(B(bidx,1)==base(S));
  t_mut(i) = sum(B(bidx,1)==base(newbase));
  others = setdiff(1:4,base([newbase S]));
  t_other1(i) = sum(B(bidx,1)==others(1));
  t_other2(i) = sum(B(bidx,1)==others(2));
end

catch me,excuse(me);end


