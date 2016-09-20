function [pred,tr,res]=read_mit_pred_results(fname)

T=read_mit_odf_file(fname);
pred=str2num(strvcat(T.data{3}));
tr=str2num(strvcat(T.data{2}));
res.conf=T.data{4};
res.correct=zeros(length(pred),1);
for i=1:length(T.data{5})
  res.correct(i)=convert_enum(T.data{5}{i},{'false',0;'true',1});
end
