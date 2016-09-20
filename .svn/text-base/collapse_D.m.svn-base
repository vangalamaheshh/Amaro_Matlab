function D=collapse_D(D,fld,method)

D=add_history(D,mfilename,fld,method);

if ischar(fld)
  if isfield(D,fld)
    f=getfield(D,fld);
    if iscell(f)
      f=strvcat(f);
    end
  else
    error('no such field');
  end
elseif size(fld,1)==size(D.dat,1)
  f=fld;
else
  f=D.gsupdat(fld,:)';
end

[u,ui,uj]=unique_keepord(f,'rows');
x=collapse_dat(D.dat,uj,1:length(ui),method);
D=reorder_D_rows(D,ui);
D.dat=x;
D.n_collapse=histc(uj,1:length(ui));
