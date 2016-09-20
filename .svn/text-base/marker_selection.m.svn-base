function [P1sgte,P1slte,P2s,Pf,rs,Pf2,S,gp,fwer,fpr,topidx,K,N]=marker_selection(dat,cls0,cls1,test_type,nperm)

D.dat=dat;
[P,S]=differential_analysis(D,cls0,cls1,test_type,0);
if nperm<0
  P2s=P;
  P1sgte=[]; P1slte=[]; Pf=[]; rs=[]; Pf2=[]; gp=[]; fwer=[]; fpr=[];
  topidx=[];
  return;
end

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
if isfield(test_type,'booster')
  nparts=1;
end

if isfield(test_type,'online')
  test_type.observed=S;
end

[PR,SR,OP,OS,rs]=differential_analysis_permutations_by_parts(D,cls0,cls1,test_type,nperm,0,nparts,lsfdir);
%Pf=fix_Pvalues(-S,-SR);
%Pf=2*min(Pf,1-Pf);
%Pf2=fix_Pvalues(-abs(S),-abs(SR));
Pf=P; %assymptotic P-value
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
%[P1sgte,P1slte,P2s]=fix_Pvalues2_by_parts(S(:,1:2),SR(:,1:2),1,0,nparts,'/xchip/data/gadgetz/lsfres/');

if isfield(test_type,'online')
  [P1sgte,P1slte,P2s]=fix_Pvalues_online(PR,SR,1);
  K=PR;
  N=SR;
else
  [P1sgte,P1slte,P2s]=fix_Pvalues2_by_parts(S,SR,1,0,nparts,['/' ...
                      'xchip/data/gadgetz/lsfres/']);
  K=[];
  N=[];
end

% save /xchip/data/gadgetz/lsfres/tmp.mat 

%[gp,fwer,fpr]=rankbased_pvalues(abs(S),abs(SR));
gp=[];
fwer=[];
fpr=[];



