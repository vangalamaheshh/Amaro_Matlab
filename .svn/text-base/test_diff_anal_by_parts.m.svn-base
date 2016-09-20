function res=test_diff_anal_by_parts(meta_dat,n,maxn,N,Nxv)

res=cell(2,Nxv);
cls0=1:(2*n);
cls1=(2*n)+(1:(2*n));
l=lsf('/xchip/data/gadgetz/lsfres/');
h=zeros(Nxv,1);
for i=1:Nxv
  D.dat=meta_dat((i-1)*N+(1:N),[1:(2*n) 2*maxn+(1:(2*n))]);
  [l,h(i)]=bsub(l,{'res_step'},'test_diff_anal_step',{D,cls0,cls1}); 
end
[l,resvec]=wait(l);
for i=1:Nxv
  res{1,i}=resvec{h(i)}.res_step{1};
  res{2,i}=resvec{h(i)}.res_step{2};
end

verbose('7')


