function [cvidx,cvq,cvs,cvr,cvpi0,Q,S]=cv_top_markers(D,supid,test_type,nperm,qthresh,frac,n)

cls0=find(D.supdat(supid,:)==0);
cls1=find(D.supdat(supid,:)==1);
n0=length(cls0);
n1=length(cls1);

nf0=round(frac*n0);
nf1=round(frac*n1);

if nf0==0 | nf1==0
  error('too low frac');
end
Q=sparse(size(D.dat,1),n);
S=Q;
for i=1:n
  r0=randperm(n0);
  r1=randperm(n1);
  cvr(i,:)=[ cls0(r0(1:nf0)) cls1(r1(1:nf1))];
  Dr=reorder_D_cols(D,cvr(i,:));
  [cvidx{i},cvq{i},cvs{i},cvpi0(i)]=get_top_markers(Dr,supid,test_type, ...
                                              nperm,qthresh);
  Q(cvidx{i},i)=cvq{i};
  S(cvidx{i},i)=cvs{i};
end

