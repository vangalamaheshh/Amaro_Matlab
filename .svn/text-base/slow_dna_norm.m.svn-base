function nv=slow_dna_norm(v,chunk)

nc=size(v,2);
for i=1:chunk:nc
  ch=i:min(i+chunk-1,nc);
  stv=std(v(ch,:),0,2);
  zero_idx=find(stv==0);		     
  if ~isempty(zero_idx)
    stv(zero_idx)=1; % if std is zero the vector is const and nv
		   % will be zero
    warning('Zero std in dna_norm');
  end
  nv(ch,:)=(v(ch,:)-repmat(mean(v(ch,:),2),1,nc))./(repmat(stv,1,nc));
end
