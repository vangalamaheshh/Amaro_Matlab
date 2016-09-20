function [nv,m,s]=dna_norm_nan(v)

nc=size(v,2);

% stv=std(v,0,2);

m=nanmean(v')';
s=nanstd(v')';
nv=(v-repmat(m,1,nc))./(repmat(s,1,nc)+eps);
if(0)
  zero_idx=find(stv==0);		     
  if ~isempty(zero_idx)
    stv(zero_idx)=1; % if std is zero the vector is const and nv
                     % will be zero
    warning('Zero std in dna_norm');
  end
  nv=(v-repmat(mean(v')',1,nc))./(repmat(stv,1,nc));
end
