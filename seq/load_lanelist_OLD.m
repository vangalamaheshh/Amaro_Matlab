function lanes = load_lanelist(fname)
if ~exist('fname','var') || isempty(fname)
  fname = '/xchip/tcga/gbm/analysis/lawrence/wgs/lanelist.txt';
end

lanes = load_struct(fname,'%f%s%s%s%s%s%f%f%s');
lanes = reorder_struct(lanes,find(~isnan(lanes.FLE)));

lanes.date = regexprep(lanes.date,'^8','08');
tmp = regexprep(lanes.TN,'T','tumor');
lanes.tumornormal = regexprep(tmp,'N','normal');

