function [D,supids]=add_sequencing_info(D,fname)

[f,fid]=read_dlm_file(fname,char(9),1);
form=repmat('%s',1,length(f{1}));

tab=read_table(fname,form,char(9),1);

nms=tab.dat{1};

[M,m1,m2]=match_string_sets(nms,D.sdesc);
if length(m1)<size(D.dat,2)
  disp(['Did not match ' num2str(size(D.dat,2)-length(m1)) ' samples']);
end

supids=[];
for i=2:3:length(tab.dat)
  v=tab.dat{i};
  v=v(m1);
  v=regexprep(v,'#N/A','NaN');
  empty_cells=find(cellfun('isempty',v));
  v(empty_cells)=cellstr(repmat('0',length(empty_cells),1));
  vn=str2num(strvcat(v));
  vn1=nan(1,size(D.dat,2));
  vn1(m2)=vn';
  [D,supid]=add_D_sup(D,tab.headers{1}{i},tab.headers{1}{i},vn1);
  supids=[supids supid];
end

