function cyto_st=pos2cyto(ch,pos,cyto)

if length(pos)>1
  st1=pos2cyto(ch,pos(1),cyto);
  st2=pos2cyto(ch,pos(end),cyto);
  if strcmp(st1,st2)
    cyto_st=st1;
  else
    cyto_st=[st1 '-' st2];
  end
  return
end
    
in_chr=find(cat(1,cyto.chrn)==ch);
cyto_st=[];
for i=1:length(in_chr)
  j=in_chr(i);
  if cyto(j).start<pos && cyto(j).end>=pos
    cyto_st=cyto(j).name;
    break;
  end
end
