function [p,vp,x,tab,msz]=sim_step(n,Fgen,Fscore,tab,ntimes,calc_pv,notab)

if ~exist('calc_pv','var')
  calc_pv=1;
end

if ~exist('notab','var')
  notab=0;
end
% make table
if (~exist('tab','var') || isempty(tab)) && ~notab
%  sz=binoinv(1-0.01/prod(size(n)),n,Fgen)+1;
%  sz=binoinv(1-1/prod(size(n)),n,Fgen)+1;
  sz=binoinv(0.99,n,Fgen)+1;
  sz=max(sz(:));
  disp(['calculating table size=' num2str(sz)]);

  if (0)
    tab=zeros([size(n) sz]);
    for i=1:sz
      tab(:,:,i)=binocdf(i-1,n,Fgen);
    end
  else
    tab=zeros([size(n) sz]);
    for i=1:sz
      tab(:,:,i)=binopdf(i-1,n,Fgen);
      if i>1
        tab(:,:,i)=tab(:,:,i-1)+tab(:,:,i);
      end
    end
  end
  if ~exist('ntimes','var')
    p=tab;
    vp=[];
    x=[];
    return
  end
end

if exist('ntimes','var') && ntimes>1
  disp(['running ' num2str(ntimes) ' times']);
  for t=1:ntimes
    if calc_pv
      [p(:,t),vp(:,t),x(:,:,t)]=sim_step(n,Fgen,Fscore,tab,1,calc_pv,notab);
    else
      [p_tmp,vp_tmp,x(:,:,t)]=sim_step(n,Fgen,Fscore,tab,1,calc_pv,notab);
      p=[];
      vp=[];
    end
    if mod(t,10)==0
      disp(t);
    end
  end
  return
end

if ~exist('Fscore','var') || isempty(Fscore)
  Fscore=Fgen;
end

r=rand(size(n));
if notab
  x=binoinv(r,n,Fgen);
else
  x=sum(repmat(r,[1 1 size(tab,3)])>tab,3);

  notintab=find(x==size(tab,3));
  if ~isempty(notintab)
    x(notintab)=binoinv(r(notintab),n(notintab),Fgen(notintab));
    warning('Reached end of table');
  end
end
%if nnz(x==size(tab,3))>0
%  error('Reached end of table');
%%   keyboard;
%end

%max(x(:))
% x=binornd(n,Fgen);
% disp('start');
if calc_pv
  T=ones(size(x));
  VT=T;
  [nzi,nzj,nzv]=find(x);
  for i=1:length(nzv)
    VT(nzi(i),nzj(i))=binopdf(nzv(i),n(nzi(i),nzj(i)),Fscore(nzi(i),nzj(i)));
    T(nzi(i),nzj(i))=1-binocdf(nzv(i)-1,n(nzi(i),nzj(i)),Fscore(nzi(i),nzj(i)));
  end
  p=prod(T,2);
  vp=prod(VT,2);
else
  p=[];
  vp=[];
end

if (0)
  T2=zeros(size(x));
  [nzi,nzj,nzv]=find(x);
  for i=1:length(nzv)
    T2(nzi(i),nzj(i))=ln_binocdf(nzv(i),n(nzi(i),nzj(i)),Fscore(nzi(i),nzj(i)),1);
  end
  lnp=sum(T,2);
end
%keyboard
% q=calc_fdr_value(p);
% sq=sort(q);


