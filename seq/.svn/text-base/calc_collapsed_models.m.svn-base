function [LL_model,refref_LL]=calc_collapsed_models(LL,refref_model,ref)

[mx_LL,mx_LLi]=max(LL,[],1);
% back to probs.
lik=10.^(LL-repmat(mx_LL,10,1)); % subtract mx_t_LL for precision issues

% AA AC AG AT CC GG TT CG CT GT
%          alt= A C G T
ref_alt_model=[ 1 2 3 4; ... % ref=A
                2 5 8 9 ; ... % ref=C  
                3 8 6 10; ... % ref=G  
                4 9 10 7; ... % ref=T  
              ];

LL_model=zeros(4,size(LL,2));
% add ref/alt
for ref_i=1:4
  idx=find(ref==ref_i);
  for alt_i=1:4
    if alt_i==ref_i, continue; end;
    LL_model(alt_i,idx)=LL_model(alt_i,idx)+lik(ref_alt_model(ref_i,alt_i),idx);
  end
end

% add alt/alt
LL_model(1:4,:)=LL_model(1:4,:)+lik(1:4,:);
LL_model=log10(LL_model)+repmat(mx_LL,4,1);

ref_N_idx=find(ref>4);
ref(ref>4)=1;

LL_model(ref+4*(0:(size(LL_model,2)-1)))=-Inf;
refref_LL=LL(refref_model+10*(0:(size(LL,2)-1)));
