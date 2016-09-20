function p_ln=hyge2pdf_ln(kv,n,k1,n1)

p_ln = nan(size(kv));
for i=1:length(kv)
    k=kv(i);
    p_ln(i)=gammaln(n1+2)-gammaln(k1+1)-gammaln(n1-k1+1) + ...
         gammaln(n+1) - gammaln(k+1) - gammaln(n-k+1) + ...
         gammaln(k1+k+1) + gammaln(n+n1-k-k1+1) - gammaln(n+n1+2);         
end

