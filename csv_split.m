function val=csv_split(txt)
if txt(1)=='"'
  txt=txt(2:end);
  end
if txt(end)=='"'
  txt=txt(1:(end-1));
end
tok=strfind(txt,'","');
tok=[-2 tok length(txt)+1];
res=[];
for i=1:(length(tok)-1)
  val{i}=txt((tok(i)+3):(tok(i+1)-1));
end
