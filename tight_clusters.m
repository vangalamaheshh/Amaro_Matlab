function [clusters,lbl]=tight_clusters(dend,level_type,level,minsz)


switch level_type
 case 'height'
  coph=lnk2cophenet(dend.lnk,0);
  G=(coph<=level);
 case 'k'
  coph=lnk2cophenet(dend.lnk,1);
  G=(coph>=k);
end  

G=triu(G,1);
[gi,gj,gk]=find(G);
lbl=unionfind([gi gj],size(G,1));
[tmp,revidx]=sort(dend.idx);
lbl=lbl(revidx);

clusters={};
for i=1:max(lbl)
  ci=find(lbl==i);
  if length(ci)>=minsz
    clusters{end+1}=ci;
  else
    break;
  end
end
