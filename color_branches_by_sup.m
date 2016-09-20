function confmat=color_branches_by_sup(Dord,supid,dend_handle,dend,nc)
coph=lnk2cophenet(dend.lnk,1);
coph(coph>=nc)=nc;
coph(coph<nc)=0;
coph=coph+eye(size(coph,1))*nc;
[dum,lbl]=max(unique(coph,'rows'),[],1);
[confmat,chi2,p,labels]=crosstab(lbl,Dord.supdat(supid,:));

disp(confmat);

empty_pos=find(cellfun('isempty',labels));
if ~isempty(empty_pos)
  labels(empty_pos)=cellstr(repmat('-1',length(empty_pos),1));
end

class_types=str2num(strvcat(labels(:,2)));
[mx,mi]=max(confmat');
color_branches(dend_handle,dend.lnk,lbl,Dord.supmark(supid).colormap(class_types(mi)+1,:));

