function D=add_gcm_scans(D,fname)
d=read_dlm_file(fname,9,2);
d1=dlmsep(d{2}{2},' ');
if isfield(D,'scans')
  D=rmfield(D,'scans');
end
D.scans=cell(size(D.sdesc,1),1);
D.scans{1}={struct('name',d1{1},'chip','HuFL'), struct('name', d1{2},'chip','Hu35KsubA')};
for i=4:3:(length(d{2})-1)
  n1=regexp(d{2}{i},'([^/]*)','match'); n1=n1{1};
  n2=regexp(d{2}{i+1},'([^/]*)','match'); n2=n2{1};
  D.scans{((i-1)/3)+1}={struct('name',n1,'chip','HuFL'), ...
                      struct('name',n2,'chip','Hu35KsubA')};
end

D.sscale=repmat('scale factor=1',size(D.dat,2),1);
