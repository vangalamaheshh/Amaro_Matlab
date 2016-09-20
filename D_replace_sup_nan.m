function D=D_replace_sup_nan(D,rc,supids)

if ~exist('supids','var') || isempty(supids)
  if is_col(rc)
    supids=1:size(D.supdat,1);
  else
    supids=1:size(D.gsupdat,1);
  end
end


if is_col(rc)
  if ischar(D.supacc)
    D.supacc=cellstr(D.supacc);
    D.supdesc=cellstr(D.supdesc);
  end
  for i=1:length(supids)
    supid=supids(i);
    if any(isnan(D.supdat(supid,:)))
      pos=find(D.supacc{supid}==':');
      if ~isempty(pos)
        n=length(find(D.supacc{supid}=='/'))+2;
        D.supdat(supid,isnan(D.supdat(supid,:)))=n;
        D.supacc{supid}=[D.supacc{supid} '/' num2str(n) '-NaN'];
        D.supdesc{supid}=[D.supdesc{supid} '/' num2str(n) '-NaN'];
      else
        D.supdat(supid,isnan(D.supdat(supid,:)))=2;
        D.supacc{supid}=[D.supacc{supid} ': 0-No/1-Yes/2-NaN'];
        D.supdesc{supid}=[D.supdesc{supid} ': 0-No/1-Yes/2-NaN'];
      end
    end
  end
else
  if ischar(D.supacc)
    D.gsupacc=cellstr(D.gsupacc);
    D.gsupdesc=cellstr(D.gsupdesc);
  end
  for i=1:length(supids)
    supid=supids(i);
    if any(isnan(D.gsupdat(supid,:)))
      pos=find(D.gsupacc{supid}==':');
      if ~isempty(pos)
        n=length(find(D.gsupacc{supid}=='/'))+2;
        D.gsupdat(supid,isnan(D.gsupdat(supid,:)))=n;
        D.gsupacc{supid}=[D.gsupacc{supid} '/' num2str(n) '-NaN'];
        D.gsupdesc{supid}=[D.gsupdesc{supid} '/' num2str(n) '-NaN'];
      else
        D.gsupdat(supid,isnan(D.gsupdat(supid,:)))=2;
        D.gsupacc{supid}=[D.gsupacc{supid} ': 0-No/1-Yes/2-NaN'];
        D.gsupdesc{supid}=[D.gsupdesc{supid} ': 0-No/1-Yes/2-NaN'];
      end
    end    
  end  
end

