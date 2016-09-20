function C=add_code_dend(C);

[ord_dat,idx]=sortrows(C.D.dat');
C.Dord=reorder_D_cols(C.D,idx);
lnk=[];
s=ones(size(ord_dat,1),1);
ord_dat=[ ord_dat (1:size(ord_dat,1))'];
for i=1:size(ord_dat,2)
  rls=runlength(s);
  for j=1:size(rls,1)
    rl=runlength(ord_dat(rls(j,1):rls(j,2),i));
    for k=1:(size(rl,1)-1)
      lnk=[lnk; (rls(j,1)-1)+rl(k,1) (rls(j,1)-1)+rl(k,2) (rls(j,1)-1)+rl(k,2)+1 rls(j,2) size(ord_dat,2)-i];
    end
    [u,ui,uj]=unique_keepord(ord_dat(:,1:i),'rows');
    s=uj';
  end
end
lnk=flipud(lnk);

C.sdend.lnk=lnk;
C.sdend.idx=idx;


if (0)
  T=C.D;
  T.dat=T.dat.*repmat(flipud(5.^(0:(size(T.dat,1)-1))'),1,size(T.dat,2));
  [Dord,dend]=one_way_clustering(T,'cols',struct('dist','euclidean','cluster','ward','preproc','none'));
  C.Dord=Dord;
  C.dend=dend;
  [sv,si]=sort(C.dend.lnk(:,5));
  C.dend.lnk(:,5)=si-1;
end
