function [supacc,supdesc]=sup_cat2name(D,supid,cat)

if ischar(supid)
  supid=find_supid(D,supid);
end

[typeacc,typedesc,D,range,non_empty]=decollapse_supdat(D,supid,1:max(D.supdat(supid,:)));
supacc=deblank(D.supacc(range(cat),:));
supdesc=deblank(D.supdesc(range(cat),:));
