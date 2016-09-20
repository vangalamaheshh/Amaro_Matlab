function svdD=svd_D(D,rc,k)

if is_col(rc)
  [u,s,v]=svd(D.dat','econ');
  uk=u(:,1:k);
  sk=s(1:k,1:k);
  vk=v(:,1:k);
  svdD=reorder_D_rows(D,[]);
  svdD.dat=uk';
  svdD.gacc=cellstr([repmat('SVD #',k,1) num2str((1:k)','%d')]);
  svdD.gdesc=svdD.gacc;
  svdD.svd.u=u;
  svdD.svd.s=s;
  svdD.svd.v=v;
else
  [u,s,v]=svd(D.dat,'econ');
  uk=u(:,1:k);
  sk=s(1:k,1:k);
  vk=v(:,1:k);
  svdD=reorder_D_cols(D,[]);
  svdD.dat=uk;
  svdD.sdesc=cellstr([repmat('SVD #',k,1) num2str((1:k)','%d')]);
  svdD.svd.m=m;
  svdD.svd.d=d;
  svdD.svd.v=v;
  svdD.svd.q=q;  
end

