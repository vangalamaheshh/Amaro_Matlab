function [Pf,PRf,Pf2]=fix_Pvalues_by_parts(P,PR,smooth_flag,nparts,lsfdir);

if nargin<3
  nparts==1;
end
if nargin<4
  lsfdir='./';
end

if nparts>1
  p=get_parts(1:size(P,1),nparts);
  l=lsf(lsfdir);
  h=zeros(nparts,1);
  for i=1:nparts
    Pi=P(p{i});
    PRi=PR(p{i},:);
    [l,h(i)]=bsub(l,{'Pf','PRf','Pf2'},'fix_Pvalues',{Pi,PRi,smooth_flag});
  end
  [l,res]=wait(l); % wait for all
  Pf=NaN*ones(size(P));
  PRf=NaN*ones(size(PR));
  Pf2=NaN*ones(size(P));
  for i=1:nparts
    Pf(p{i})=res{h(i)}.Pf;
    PRf(p{i},:)=res{h(i)}.PRf;
    Pf2(p{i},:)=res{h(i)}.Pf2;
  end
else
  [Pf,PRf,Pf2]=fix_Pvalues(P,PR,smooth_flag);
end

