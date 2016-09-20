function C=determine_euploid_level_and_stromal_contamination(C,mtype);

if ischar(mtype)
  tmp.method=mtype;
  mtype=tmp;
end

for i=1:size(C.dat,2)
  switch mtype.method
   case 'fit'
    if isfield(C,'peaks_allele_balance')
      ab=C.peaks_allele_balance{i};
    else
      ab=ones(1,length(C.peaks{i}));
    end
    hv=histc(C.level(:,i),1:length(C.peaks{i}));
    hv=as_row(hv)./sum(hv);
    
    cand=find(hv>mtype.min_frac & ab>mtype.min_allele_balance);
    C.trans{i}.CAND=cand;
    
    %FIXME: not should be like this.
    if isempty(cand)
       [tmp,euploid]=max(hv);
    else
      [tmp,euploid]=max(hv(cand));
      euploid=cand(euploid);
    end
    
    C.trans{i}.euploid=euploid;
    if isempty(cand)
      [tmp,cand]=max(hv);
    end
    
    for j=1:length(C.peaks{i})
      in_pk=find(C.level(:,i)==j);
      C.peaks_se{i}(j)=mad(C.raw(in_pk,i))*mad_factor;
    end
    
    r=2.^(C.peaks{i}+1);
    dr=2.^(C.peaks{i}+1).*C.peaks_se{i}*log(2);
    for ei=1:length(mtype.euploid_levels)
      eu=mtype.euploid_levels(ei);
      for ci=1:length(cand)
        [p,alpha,s2,cn,pval]=stromal_cont_estimation(r,cand(ci),eu,[-1 1]);
        C.trans{i}.S2(ci,ei)=s2;
        C.trans{i}.P(ci,ei)=p;
        C.trans{i}.PVAL(ci,ei,:)=pval;
        C.trans{i}.CN(ci,ei,:)=cn;
      end
    end
    
    if (0)
      g=zeros(length(mtype.purity_levels),length(cand),length(mtype.euploid_levels));
      for ei=1:length(mtype.euploid_levels)
        eu=mtype.euploid_levels(ei);
        for ci=1:length(cand)
          alpha=eu
          for pi=1:length(mtype.purity_levels)
            p=mtype.purity_levels(pi);
            
          end
        end
      end
    end
    
  end
end
