function D=light_D(D,dirname,fname,ext,flds,v)

for i=1:length(flds)
  if isfield(D,flds{i})
    if ~isfield(D,'light_flds')
      D.light_flds=[];
      last_fld=0;
    else
      last_fld=length(D.light_flds);
    end
    D.light_flds(last_fld+1).name=flds{i};
    x=getfield(D,flds{i});
    outname=[ dirname fname ext '.' flds{i} '.mat'];
    if exist('v','var')
      save(outname,'x',v);
    else
      save(outname,'x');
    end
    D.light_flds(last_fld+1).fname=outname;
    D=setfield(D,flds{i},[]);
  end
end
