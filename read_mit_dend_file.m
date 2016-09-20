function dend=read_mit_dend_file(fname)

odf=read_mit_odf_file(fname);
dend.Model=odf.Model;
dend.idx=str2num(odf.NewOrder);
dend.lnk=zeros(length(odf.data{1}),5);
for i=1:5
  dend.lnk(:,i)=double(odf.data{i});
end
[dend.g,dend.v,dend.m]=lnk_to_graph(dend.lnk,length(dend.idx));
[dum,dend.revidx]=sort(dend.idx);
dend.z=lnk2z(dend.lnk,dend.idx);

