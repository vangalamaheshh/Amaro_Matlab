function write_mit_dend_file(fname,dend)

odf.NewOrder=sprintf('%d ',dend.idx);
odf.col_names={'Left from','Left to','Right from','Right to', ...
               'Value'};
odf.col_types={'int','int','int','int','float'};
for i=1:5
  odf.data{i}=dend.lnk(:,i);
end
if isfield(dend,'Model')
  odf.Model=dend.Model
else
  odf.Model='Dendrogram';
end

write_mit_odf_file(fname,odf);

