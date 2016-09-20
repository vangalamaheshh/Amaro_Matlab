function lhs=show_dend_lnk_idx(lnk,idx,d,sz,gtype,sc,lw)

if ~exist('lw','var')
  lw=[];
end

N=length(idx);

if ~isempty(d)  
  [idx2,lnk2]=reorder_dend(d(idx,idx),lnk);
  idx2=idx(idx2);
else
  idx2=idx;
  lnk2=lnk;
end


[g,v,m]=lnk_to_graph(lnk2,N);

if lnk(1,5)<lnk(end,5)
  nan_val=0;
else
  nan_val=max(lnk2(:,5)+1); % for SPC
end

s=ones(size(m,1),1);
if ~isempty(sc)
  s(1:N)=sc(idx2);
  s((N+1):end)=-1;
end

small=find(m(:,2)-m(:,1)+1 < sz );
big=find(m(:,2)-m(:,1)+1 >= sz );
all_lnk=1:size(m,1);
g1=g;
g1(small,:)=0;
g2=reduce_graph(g1);
% keyboard
if (nargin <= 7) | (noshow==0)
  switch (gtype)
   case {'all'}
    lhs=show_dend_graph(g,v,m,s,N,nan_val,lw);
   case {'big'}
    lhs=show_dend_graph(g1,v,m,s,N,nan_val,lw);  
   case {'reduced_big'}
    lhs=show_dend_graph(g2,v,m,s,N,nan_val,lw);  
  end
end

