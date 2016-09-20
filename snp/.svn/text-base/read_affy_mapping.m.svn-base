function res=read_affy_mapping(fname)

res=read_table(fname,[],',',1,'emptyvalue',NaN,'bufsize',1e6);

[ln,fid]=read_dlm_file(fname,[],1);
res.headers{1}=csv_split(ln{1}{1});

tmp=textscan(fid,'%s','delimiter','\n','bufsize',1e6);
tmp=tmp{1};
for i=1:length(tmp)
  res.dat{i}=csv_split(tmp{i});
end
res.dat=cat(1,res.dat{:});

