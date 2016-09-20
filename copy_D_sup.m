function [Dto,r]=copy_D_sup(Dfrom,Dto,range,rc)

if ~exist('rc','var')
  rc='col';
end

if ~exist('range','var') || ( exist('range','var') && ...
                              isempty(range))
  if is_col(rc)
    range=size(Dfrom.supdat,1);
  else
    range=size(Dfrom.gsupdat,1);
  end
end

if is_col(rc)
  Dto.supdat=[Dto.supdat; Dfrom.supdat(range,:)];
  if iscell(Dto.supacc)
    Dto.supacc=[Dto.supacc; Dfrom.supacc(range)];
  else
    Dto.supacc=strvcat(Dto.supacc,deblank(Dfrom.supacc(range,:)));
  end
  if iscell(Dto.supdesc)
    Dto.supdesc=[Dto.supdesc; Dfrom.supdesc(range)];
  else
    Dto.supdesc=strvcat(Dto.supdesc,deblank(Dfrom.supdesc(range,:)));
  end
  r=(size(Dto.supdat,1)-length(range)+1):size(Dto.supdat,1);
else
  Dto.gsupdat=[Dto.gsupdat; Dfrom.gsupdat(range,:)];
  if iscell(Dto.gsupacc)
    Dto.gsupacc=[Dto.gsupacc; Dfrom.gsupacc(range)];
  else
    Dto.gsupacc=strvcat(Dto.gsupacc,deblank(Dfrom.gsupacc(range,:)));
  end
  if iscell(Dto.gsupdesc)
    Dto.gsupdesc=[Dto.gsupdesc; Dfrom.gsupdesc(range)];
  else
    Dto.gsupdesc=strvcat(Dto.gsupdesc,deblank(Dfrom.gsupdesc(range,:)));
  end
  r=(size(Dto.supdat,1)-length(range)+1):size(Dto.supdat,1);
end  
