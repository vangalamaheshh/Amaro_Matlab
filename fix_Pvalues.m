function [Pf,PRf,Pf2]=fix_Pvalues(P,PR,smooth_flag);

if nargin==2
  smooth_flag=0;
end

Pf=P;
Pf2=P;
PRf=PR;
n=size(PR,2);
[sv,si]=sort(PR');
[ssi,ri]=sort(si);

for i=1:length(Pf)
%  [sv,si]=sort(PR(i,:));
%  [ssi,ri]=sort(si);
  
  tr_less=length(find(sv(:,i)<P(i)));
  tr_equ=length(find(sv(:,i)==P(i)));
  tr_more=length(find(sv(:,i)>P(i)));
  
  cn=tr_less+tr_equ+tr_more; % ignore NaNs
  if tr_less<tr_more
    tr1=tr_less;
    tr=tr_less+tr_equ;
    tr2=tr_more;
    other_tail_equ=length(find(sv((tr+1):end,i)==sv(end-tr+1,i)));
    other_tail_more=length(find(sv((tr+1):end,i)>sv(end-tr+1,i)));
    if (other_tail_equ+other_tail_more)>tr
      other_tail=other_tail_more;
    else
      other_tail=other_tail_equ+other_tail_more;
    end
  else
    tr1=tr_more;
    tr=tr_more+tr_equ;
    tr2=tr_less;
    other_tail_equ=length(find(sv(1:(tr-1),i)==sv(tr,i)));
    other_tail_less=length(find(sv(1:(tr-1),i)<sv(tr,i)));
    if (other_tail_equ+other_tail_less)>tr
      other_tail=other_tail_less;
    else
      other_tail=other_tail_equ+other_tail_less;
    end
  end
  if (smooth_flag)
%    Pf(i)=(tr+1)/(cn+1);
%    PRf(i,:)=(ri(:,i)+1)/(cn+1);
    Pf(i)=(tr+1)/(cn+2);
    if (tr>0) & (tr<size(sv,1))
      Pf2(i)=(tr+other_tail+1)/(cn+2);
    else
      Pf2(i)=Pf(i);
    end
    PRf(i,:)=(ri(:,i)+1)/(cn+2);
  else
    Pf(i)=tr/cn;
    Pf2(i)=(tr+other_tail)/cn;
    PRf(i,:)=ri(:,i)/cn;
  end
  if mod(i,1000)==0
    verbose(i);
  end
end
