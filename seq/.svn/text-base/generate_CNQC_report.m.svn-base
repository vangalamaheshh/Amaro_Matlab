function L = generate_CNQC_report(Lt,Ln,individual_name)
% Mike Lawrence 2009-10-28

if ~exist('individual_name','var'), individual_name=''; end

L = [];
L.individual = repmat({individual_name},slength(Lt)+slength(Ln),1);
L.tn = [repmat({'tumor'},slength(Lt),1);repmat({'normal'},slength(Ln),1)];
flds1 = {'SM','ID','PU','CN','PL','DT','LB','baitset',...
  'cov','enoughreads','noise','normalness','tumorseg_corr','tumormed_corr','judgement','is_mixup','is_blacklisted'};
flds2 = {'sample','readgroup','flowcell_lane','center','platform','date','library','baitset',...
  'nreads','enoughreads','noise','normalness','tumorseg_corr','tumormed_corr','judgement','is_mixup','is_blacklisted'};
if ~isfield(Lt,'ID'), Lt.ID = repmat({'-'},slength(Lt),1); end
if ~isfield(Lt,'baitset'), Lt.baitset = repmat({'-'},slength(Lt),1); end
if ~isfield(Ln,'ID'), Ln.ID = repmat({'-'},slength(Ln),1); end
if ~isfield(Ln,'baitset'), Ln.baitset = repmat({'-'},slength(Ln),1); end
tmp = concat_structs({keep_fields(Lt,flds1),keep_fields(Ln,flds1)});
tmp = keep_fields(tmp,flds1);
tmp = rename_fields(tmp,flds1,flds2);
L = merge_structs({L,tmp});
