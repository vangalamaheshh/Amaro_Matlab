function C=make_cipher(D,use_base_entropy)

if ~exist('use_base_entropy','var')
  use_base_entropy=1;
end

if use_base_entropy
  disp('using base_entropy');
end

C=[];
N=size(D.dat,2);

c=zeros(1,N);
h(1)=calc_entropy(c);

% n0*(n1+n2)+n1*n2 = n0*n1+n0*n2+n1*n2
% n0*log(n0)+n1*log(n1)+n2*log(n2)
D.gsupacc='BE';
D.gsupacc='BaseEnt';
for i=1:size(D.dat,1)
  D.gsupdat(1,i)=calc_entropy(D.dat(i,:));
end

Dsave=D;

D.gorigidx=1:size(D.dat,1);


tc=c;
gorigidx=[];
while (size(D.dat,1)>0) && (h(end)>0) && ((length(h)==1) ||  (h(end)~=h(end-1))) 
  ng=size(D.dat,1);
  th=zeros(ng,1);
  for i=1:ng;
    tc=[c; D.dat(i,:)];
    th(i)=calc_entropy(tc);
  end
  minh=min(th);
  idx=find(th==minh);
  if length(idx)>1
    disp(['Chosing min base entropy among ' num2str(length(idx))]);
    if use_base_entropy
      [mn,mni]=min(D.gsupdat(1,idx),[],2);
      minbei=find(D.gsupdat(1,idx)==mn);
      if length(minbei)>1
        disp(['Chosing min BE among ' num2str(length(minbei))]);
      end
      if (mni~=1)
        disp(['Delta BE: ' num2str(mn-D.gsupdat(1,idx(1)))]);
      end
      idx=idx(mni);
    else
      idx=idx(1);
    end
  end
  c=[ c; D.dat(idx,:)];
  gorigidx=[gorigidx; D.gorigidx(idx)];
  h=[h; th(idx)];
  D=reorder_D_rows(D,setdiff(1:size(D.dat,1),idx));
end
C.h=h;
C.code=c(2:end,:);
C.idx=gorigidx;
C.D=rmfield(reorder_D_rows(Dsave,C.idx),'orig');


function h=calc_entropy(code)
[u,ui,uj]=unique(code','rows');
hc=histc(uj,1:max(uj));
h=1/size(code,2)*sum(hc.*log2(hc),1);

