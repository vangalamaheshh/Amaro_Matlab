function p=hyge2pdf_at_zero(n,k1,n1)

%p = nan(size(kv));
%for i=1:length(kv)
%    k=kv(i);
%    term(i)=gammaln(n1+2)-gammaln(k1+1)-gammaln(n1-k1+1) + ...
%         gammaln(n+1) - gammaln(k+1) - gammaln(n-k+1) + ...
%         gammaln(k1+k+1) + gammaln(n+n1-k-k1+1) - gammaln(n+n1+2);
%    p(i)=exp(term(i));
%end
%if all(kv-(0:n)==0)
%    p=p./sum(p);
%end

%term = gammaln(n1+2)-gammaln(k1+1)-gammaln(n1-k1+1) + ...
%        gammaln(n+1) - gammaln(1) - gammaln(n+1) + ...
%        gammaln(k1+1) + gammaln(n+n1-k1+1) - gammaln(n+n1+2);         

method = 1;

if method==1

  term = gammaln(n1+2)-gammaln(n1-k1+1) + ...
       + gammaln(n+n1-k1+1) - gammaln(n+n1+2);

  p = exp(term);

elseif method==2

  % METHOD 2 DOESN'T WORK YET

  if length(n)>1 || length(k1)>1 || length(n1)>1
    error('method_2 not vectorized');
  end

  rr = geoseries(1e-100,0.99,1000);
  lb = betapdf(rr,k1+1,n1+1);
  lbnorm = lb/sum(lb);
  prob_bin = 1-((1-rr).^n);
  post = lbnorm.*prob_bin;

  pr(rr,lb,lbnorm,prob_bin,post,cumsum(post));

  p = sum(post);

end


p(p<0) = 0;
p(p>1) = 1;
