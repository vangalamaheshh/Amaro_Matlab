function vec=ind2sub_vec(sz,ind)
idx=repmat(ind,1,length(sz));
cp=cumprod([1 sz]);
cp=cp(1:end-1);
q1=floor((idx-1)./cp);
q1=q1-floor(q1./sz).*sz;
vec=q1+1;
