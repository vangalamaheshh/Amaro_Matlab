function [P1sgte,P1slte,P2s]=fix_Pvalues2_by_parts(S,SR,smooth_flag,trueperm_is_first,nparts,lsfdir);

if nargin<3
  nparts==1;
end
if nargin<4
  lsfdir='./';
end

if nparts>1
  p=get_parts(1:size(S,1),nparts);
  l=lsf(lsfdir);
  h=zeros(nparts,1);
  for i=1:nparts
    Si=S(p{i});
    SRi=SR(p{i},:);
    [l,h(i)]=bsub(l,{'P1sgte','P1slte','P2s'},'fix_Pvalues2',{Si,SRi,smooth_flag,trueperm_is_first});
  end
  szSR=size(SR);
  clear SR;
  [l,res]=wait(l); % wait for all
  clear l;
  P1sgte=NaN*ones(szSR+(1-trueperm_is_first)*[0 1]);
  P1slte=NaN*ones(szSR+(1-trueperm_is_first)*[0 1]);
  P2s=NaN*ones(szSR+(1-trueperm_is_first)*[0 1]);
  for i=1:nparts
    P1sgte(p{i},:)=res{h(i)}.P1sgte;
    P1slte(p{i},:)=res{h(i)}.P1slte;
    P2s(p{i},:)=res{h(i)}.P2s;
    red{h(i)}=[];
  end
else
  [P1sgte,P1slte,P2s]=fix_Pvalues2(S,SR,smooth_flag,trueperm_is_first);
end

