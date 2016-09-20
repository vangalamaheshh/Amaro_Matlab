function [s,D1]=get_D_scans(D,chiptype)
% get the scans from type chiptype

s=[];
notempty=find(~cellfun('isempty',D.scans));
included=[];
for i=1:length(notempty)
  i1=notempty(i);
  for j=1:length(D.scans{i1})
    if strcmp(deblank(D.scans{i1}{j}.chip),chiptype)
      s=strvcat(s,D.scans{i1}{j}.name);
      included=[included i1];
    else
      verbose(['* Did not include ' D.scans{i1}{j}.chip ' ' D.scans{i1}{j}.name]);
    end
  end
end
D1=reorder_D_cols(D,included);


