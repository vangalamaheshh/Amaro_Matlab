function [scale_factor,invarient_set]=invarient_set_normalization(D,threshold)

D1=filter_D_rows(D,struct('method','minval','thresh',32));
R=D1;
for i=1:max(D.gsupdat(2,:))
  bsi=find(D1.gsupdat(2,:)==i);
%   disp([ i length(bsi)]);
  [tmp,ri]=sort(D1.dat(bsi,:),1);
  [tmp,r]=sort(ri,1);
  R.dat(bsi,:)=r;
end

s=std(R.dat,0,2);
low_s=find(s<threshold);

[Mt,m1,m2]=match_string_sets(R.gacc(low_s),D.gacc);

x=log2(D.dat(m2,:));
for i=1:max(D.gsupdat(2,:))
  bs_i=find(R.gsupdat(2,low_s)==i);
  if isempty(bs_i)
    error('no beads for this beadset');
  else
    tot(i,:)=mean(x(bs_i,:),1);
  end
end

tmp=median(sum(tot),2);
[tmp2,scale_to]=min(abs(sum(tot)-tmp));

scale_factor=2.^(repmat(tot(:,scale_to),1,size(tot,2))-tot);
invarient_set=m2;
