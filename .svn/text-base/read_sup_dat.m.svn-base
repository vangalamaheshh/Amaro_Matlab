function [supdat,sup_desc,sup_acc,cdesc]=read_sup_dat(fname)
l=my_textread(fname,'%s','whitespace','\n');
remove=[];
for i=1:length(l)
  if l{i}(1)=='#'  % # in first char. is a comment line
    remove=[remove i];
  end
end
l=l(setdiff(1:length(l),remove));

pos=findstr(9,l{1});
nc=size(pos,2);
cdesc=[];
for j=3:nc
   cdesc=strvcat(cdesc, l{1}((pos(j-1)+1):(pos(j)-1)));
end
cdesc=strvcat(cdesc, l{1}((pos(j)+1):end));

sup_desc=[];
sup_acc=[];
nr=size(l,1)-2;
supdat=zeros(nr,nc-1);
for k=2:size(l,1)
    pos=findstr(9,l{k});
    sup_acc=strvcat(sup_acc,l{k}(1:pos(1)-1));
    sup_desc=strvcat(sup_desc,l{k}((pos(1)+1):(pos(2)-1)));
    for j=2:nc
      if (j<nc)
	s=l{k}((pos(j)+1):(pos(j+1)-1));
      else
	s=l{k}((pos(j)+1):end);
      end
      if isempty(s)
         supdat(k-1,j-1)=NaN;
      else
         supdat(k-1,j-1)=sscanf(s,'%f');
      end
   end      
end





