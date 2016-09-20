function p=hyge2pdf(k,n,k1,n1)

%p = nan(size(kv));
%for i=1:length(kv)
%    k=kv(i);
%    term(i)=gammaln(n1+2)-gammaln(k1+1)-gammaln(n1-k1+1) + ...
%         gammaln(n+1) - gammaln(k+1) - gammaln(n-k+1) + ...
%         gammaln(k1+k+1) + gammaln(n+n1-k-k1+1) - gammaln(n+n1+2);         
%    p(i)=exp(term(i));
%end

term = gammaln(n1+2)-gammaln(k1+1)-gammaln(n1-k1+1) + ...
        gammaln(n+1) - gammaln(k+1) - gammaln(n-k+1) + ...
        gammaln(k1+k+1) + gammaln(n+n1-k-k1+1) - gammaln(n+n1+2);
p = exp(term);

p(p<0) = 0;
p(p>1) = 1;

%if all(k-(0:n)==0)
%    p=p./sum(p);
%end
