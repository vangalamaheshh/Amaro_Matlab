function stc=splitstr(str,sep)

spos=findstr(str,sep);
if isempty(spos)
   stc={str};
else
   spos=[spos length(str)+1];
   stc=cell(length(spos),1);
   p=1;
   for i=1:length(spos)
      stc{i}=str(p:(spos(i)-1));
      p=spos(i)+length(sep);
   end
end
