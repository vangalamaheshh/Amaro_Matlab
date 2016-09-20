function [idx,q,p,s,pi0,F]=get_top_markers(D,supid,test_type,nperm,selection)

if exist('my_gp','var')
  [r,fet,s,p,fpr,fwer,rankp,fdr,q,pi0,T]=gp_marker_selection(D,supid, ...
                                                    test_type, ...
                                                    nperm);
  p2=(p*nperm+1)/(nperm+1);
  p2(find(isnan(s)))=NaN;
else
  u=unique(D.supdat(supid,:)); 
  u(isnan(u))=[];
  assert(length(u)==2,'should be only 2 types');
  [P1sgte,P1slte,P2s,Pf,rs,Pf2,S,gp,fwer,fpr]=marker_selection(D.dat,...
                                                    find(D.supdat(supid,:)==u(1)),...
                                                    find(D.supdat(supid,:)==u(2)),...
                                                    test_type, ...
                                                    nperm);
  if size(P2s,2)>1
    p2=P2s(:,1);
    s=S(:,1);
  else
    p2=P2s;
    s=S;
  end
end
[q2,pi02]=calc_q_value(p2);
pi0=pi02;
fdr2=calc_fdr_value(p2);
switch selection.method
 case 'top'
  [sp2,sp2i]=sort(p2);
  idx=sp2i(1:min(selection.n,length(p2)));
 case 'bhfdr'
  idx=find(fdr2<=selection.thresh);
 case 'stfdr'
  idx=find(q2<=selection.thresh);
 case 'bonferroni'
  idx=find(p2*length(p2)<=selection.thresh);
 otherwise
  error('no such method');
end
F.q=q2;
F.p=p2;
F.s=s;
F.pi0=pi02;
F.fdr=fdr2;
if ~isempty(idx)
  q=q2(idx);
  s=s(idx);
  p=p2(idx);
else
  s=[];
  q=[];
  p=[];
end

