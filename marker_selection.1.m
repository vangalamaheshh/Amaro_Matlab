function [P1sgte,P1slte,P2s,Pf,rs,Pf2,S,gp,fwer,fpr,topidx]=marker_selection(dat,cls0,cls1,test_type,nperm)

D.dat=dat;
[P,S]=differential_analysis(D,cls0,cls1,test_type,1);
if ~ischar(test_type) && isfield(test_type,'nparts_perm')
  nparts=test_type.nparts_perm;
else
  nparts=10;
end
if isstruct(test_type) && isfield(test_type,'calc_p_for_top')
  [dum,topidx]=sort(abs(S));
  topidx=topidx((end-test_type.calc_p_for_top+1):end);
  D.dat=D.dat(topidx,:); % should use reorder_D_rows(D,topidx);
else
  topidx=1:size(D.dat,1);
end
S=S(topidx);

lsfdir='/xchip/data/gadgetz/lsfres/';
[PR,SR,OP,OS,rs]=differential_analysis_permutations_by_parts(D,cls0,cls1,test_type,nperm,0,nparts,lsfdir);
%Pf=fix_Pvalues(-S,-SR);
%Pf=2*min(Pf,1-Pf);
%Pf2=fix_Pvalues(-abs(S),-abs(SR));
Pf=[];
Pf2=[];

if prod(size(SR))>1e5
  if ~ischar(test_type) && isfield(test_type,'nparts_fix')
    nparts=test_type.nparts_fix;
  else
    nparts=10;
  end
else
  nparts=1;
end
[P1sgte,P1slte,P2s]=fix_Pvalues2_by_parts(S,SR,1,0,nparts,'/xchip/data/gadgetz/lsfres/');

save /xchip/data/gadgetz/lsfres/tmp.mat 

%[gp,fwer,fpr]=rankbased_pvalues(abs(S),abs(SR));
gp=[];
fwer=[];
fpr=[];



