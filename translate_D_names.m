function D=translate_D_names(D,old_names,new_names,rc)

if ~exist('rc','var')
  rc='col';
end

switch rc(1:min(length(rc),3))
 case {'col','con','exp','sam','sde'}
  [idx,had_a_match]=get_D_idx(D,'col',old_names);
  nms=cellstr(D.sdesc);
  nms(idx)=new_names(had_a_match);
  D.sdesc=strvcat(nms);
 case {'row','gen','mir','gac'}
  [idx,had_a_match]=get_D_idx(D,'row',old_names);
  if ~iscell(D.gacc)
    nms=cellstr(D.gacc);
  else
    nms=D.gacc;
  end
  nms(idx)=new_names(had_a_match);
  D.acc=nms;
 case {'gde'}
  [idx,had_a_match]=get_D_idx(D,'row',old_names);
  if ~iscell(D.gdesc)
    nms=cellstr(D.gdesc);
  else
    nms=D.gdesc;
  end
  nms(idx)=new_names(had_a_match);
  D.gdesc=nms;  
end

  
