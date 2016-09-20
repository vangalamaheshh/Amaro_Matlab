function pcaD=pca_D(D,rc,k)

if is_col(rc)
  [prj,m,d,v,q]=pca(D.dat',k);
  pcaD=reorder_D_rows(D,[]);
  pcaD.dat=prj';
  pcaD.gacc=cellstr([repmat('PC #',k,1) num2str((1:k)','%d')]);
  pcaD.gdesc=pcaD.gacc;
  pcaD.pca.m=m;
  pcaD.pca.d=d;
  pcaD.pca.v=v;
  pcaD.pca.q=q;  
else
  [prj,m,D,V,Q]=pca(D.dat,k);
  pcaD=reorder_D_cols(D,[]);
  pcaD.dat=prj;
  pcaD.sdesc=cellstr([repmat('PC #',k,1) num2str((1:k)','%d')]);
  pcaD.pca.m=m;
  pcaD.pca.d=d;
  pcaD.pca.v=v;
  pcaD.pca.q=q;  
end
